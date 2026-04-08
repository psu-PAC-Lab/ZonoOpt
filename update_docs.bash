# make C++ documentation
cd docs/C++
git rm -rf html
doxygen
git add html
cd ..

# make python documentation
cd python
git rm -rf build
make html
sed -i 's/images/_static/g' build/html/index.html
git add build
cd ../..