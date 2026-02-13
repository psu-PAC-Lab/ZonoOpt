#ifndef ZONOOPT_MI_DATA_STRUCTURES_
#define ZONOOPT_MI_DATA_STRUCTURES_

/**
 * @file MI_DataStructures.hpp
 * @author Josh Robbins (jrobbins@psu.edu)
 * @brief Data structures for mixed-integer optimization in ZonoOpt library.
 * @version 1.0
 * @date 2025-06-04
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <queue>
#include <mutex>
#include <utility>
#include <memory>
#include <set>
#include <algorithm>
#include <vector>
#include <functional>

#include "ADMM.hpp"
#include "SolverDataStructures.hpp"

namespace ZonoOpt::detail
{
    struct MI_data
    {
        std::shared_ptr<ADMM_data> admm_data;
        std::pair<int, int> idx_b; // indices of binary variables
        bool zero_one_form = false; // true: binary->{0,1}. false: binary->{-1,1}
    };

    template <typename T>
    class ThreadSafeAccess
    {
    public:
        ThreadSafeAccess() = default;

        // get method
        T get() const
        {
            std::lock_guard<std::mutex> lock(mtx);
            return data;
        }

        // set method
        void set(const T& value)
        {
            std::lock_guard<std::mutex> lock(mtx);
            data = value;
        }

    private:
        mutable std::mutex mtx;
        T data;
    };

    template <typename T>
    class ThreadSafeIncrementable
    {
    public:
        explicit ThreadSafeIncrementable(T value) : data(value)
        {
        }

        // get method
        T get() const
        {
            std::lock_guard<std::mutex> lock(mtx);
            return data;
        }

        // set method
        void set(T value)
        {
            std::lock_guard<std::mutex> lock(mtx);
            data = value;
        }

        // increment
        void operator+=(T value)
        {
            std::lock_guard<std::mutex> lock(mtx);
            data += value;
        }

    private:
        mutable std::mutex mtx;
        T data;
    };

    class ThreadSafeMultiset
    {
    public:
        // constructor
        ThreadSafeMultiset() = default;

        // add J
        void add(const zono_float J)
        {
            std::lock_guard<std::mutex> lock(mtx);
            data.insert(J);
        }

        // remove J
        void remove(const zono_float J)
        {
            std::lock_guard<std::mutex> lock(mtx);
            if (const auto it = data.find(J); it != data.end())
                data.erase(it);
        }

        // clear
        void clear()
        {
            std::lock_guard<std::mutex> lock(mtx);
            data.clear();
        }

        // get min, bool indicates if valid (i.e. not empty)
        std::pair<zono_float, bool> get_min() const
        {
            std::lock_guard<std::mutex> lock(mtx);
            if (data.empty())
                return std::make_pair(std::numeric_limits<zono_float>::infinity(), false);
            else
                return std::make_pair(*data.begin(), true);
        }

        // get size
        size_t size() const
        {
            std::lock_guard<std::mutex> lock(mtx);
            return data.size();
        }

    private:
        mutable std::mutex mtx;
        std::multiset<zono_float> data;
    };

    template <typename T>
    class ThreadSafeVector
    {
    public:
        // constructor
        ThreadSafeVector() = default;

        // add element
        void push_back(const T& value)
        {
            std::lock_guard<std::mutex> lock(mtx);
            data.push_back(value);
        }

        // clear vector
        void clear()
        {
            std::lock_guard<std::mutex> lock(mtx);
            data.clear();
        }

        // get size
        size_t size() const
        {
            std::lock_guard<std::mutex> lock(mtx);
            return data.size();
        }

        // get vector
        std::vector<T> get() const
        {
            std::lock_guard<std::mutex> lock(mtx);
            return data;
        }

        // check containment
        bool contains(const T& value, std::function<bool(const OptSolution&, const OptSolution&)>& compare) const
        {
            std::lock_guard<std::mutex> lock(mtx);
            for (auto it = data.begin(); it != data.end(); ++it)
            {
                if (compare(*it, value)) return true;
            }
            return false;
        }

    private:
        mutable std::mutex mtx;
        std::vector<T> data;
    };

    class Node final : public ADMM_solver
    {
    public:
        // constructors
        explicit Node(const std::shared_ptr<ADMM_data>& data) : ADMM_solver(data)
        {
            // init box constraints
            this->x_box = *data->x_box;
        }

        // copy constructor
        Node(const Node& other) : ADMM_solver(other)
        {
            this->x_box = other.x_box;
            this->solution = other.solution;
            this->width = this->x_box.width();
        }

        // public fields
        OptSolution solution = OptSolution();
        zono_float width = std::numeric_limits<zono_float>::infinity(); // width of the box

        // run contractor
        bool run_contractor()
        {
            const bool contractor_feasible = this->startup(this->x_box, this->solution);
            this->width = this->x_box.width();
            return contractor_feasible;
        }

        // fix bound
        bool fix_bound(int ind, const zono_float val)
        {
            // update bounds
            this->x_box[ind] = Interval(val, val);

            // run contractor
            const bool contractor_feasible = this->startup(this->x_box, this->solution, {ind});
            this->width = this->x_box.width();
            return contractor_feasible;
        }

        // solve
        void solve(std::atomic<bool>* stop = nullptr)
        {
            this->solve_core(this->x_box, this->solution, stop);
        }

        // update convergence tolerances
        void update_convergence_tolerances(const zono_float eps_prim, const zono_float eps_dual)
        {
            this->eps_prim = eps_prim;
            this->eps_dual = eps_dual;
        }

        // get box
        const Box& get_box() const
        {
            return x_box;
        }

    private:
        // members
        Box x_box;
    };

    template <typename T, typename Compare=std::less<T>>
    class PriorityQueuePrunable final : public std::priority_queue<T, std::vector<T>, Compare>
    {
    public:
        // constructor
        explicit PriorityQueuePrunable(const Compare& comp = Compare()) : std::priority_queue<
            T, std::vector<T>, Compare>(comp)
        {
        }

        // prune queue
        void prune(const T& t)
        {
            auto it_prune = std::find_if(this->c.begin(), this->c.end(), [&](const T& item)
            {
                return this->comp(item, t);
            });
            if (it_prune != this->c.end())
            {
                this->c.erase(it_prune, this->c.end());
            }
        }

        // bottom
        const T& bottom() const
        {
            return this->c.back();
        }

        // pop top
        T pop_top()
        {
            T top = std::move(this->c.front());
            this->pop();
            return top;
        }

        // clear
        void clear()
        {
            while (!this->empty())
                this->pop();
        }
    };
} // namespace ZonoOpt::detail


#endif
