rm -rf build

swig -c++ -python Aono.i #1> /dev/null 2> /dev/null
python setup.py install #1> /dev/null 2> /dev/null
