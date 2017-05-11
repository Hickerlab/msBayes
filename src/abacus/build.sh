#!/bin/sh

usage () {
    echo ""
    echo "usage: $0 [-h|--help] [-p|--prefix <INSTALL-PREFIX>] [-s|--static] [-i|--install] [-t|--test]"
    echo "  -h|--help       show help message and exit."
    echo "  -s|--static     build statically linked binaries."
    echo "  -p|--prefix     install path. Default: '/usr/local'."
    echo "  -t|--test       run the test suite after building."
    echo "  -i|--install    install the executables."
    echo ""
}

# get location of script
this_dir=`dirname "$0"`
if [ "$this_dir" = "." ]
then
    base_dir="$(pwd)"
else
    cd "$this_dir"
    base_dir="$(pwd)"
fi

# process args
extra_args=""
static=""
install_prefix=""
run_tests=""
run_install=""
while [ "$1" != "" ]
do
    case $1 in
        -h| --help)
            usage
            exit
            ;;
        -p| --prefix)
            shift
            install_prefix=$1
            ;;
        -s| --static)
            static=1
            ;;
        -t| --test)
            run_tests=1
            ;;
        -i| --install)
            run_install=1
            ;;
        * )
            extra_args="$extra_args $1"
    esac
    shift
done
args=""
if [ -n "$install_prefix" ]
then
    args="${args} -DCMAKE_INSTALL_PREFIX=${install_prefix}"
fi
if [ -n "$static" ]
then
    args="${args} -DSTATIC_LINKING=ON"
fi
if [ -n "$extra_args" ]
then
    args="${args} ${extra_args}"
fi

# check for build directory
build_dir="${base_dir}/build"
if [ -d "$build_dir" ]
then
    echo "ERROR: build directory '$build_dir' already exists."
    echo "To reconfigure, please remove this directory and re-run this script."
    exit 1
else
    mkdir "$build_dir"
fi

# configure make files and build
cd "$build_dir"
echo "Configuring make files..."
echo "    cmake ../ $args"
cmake ../ $args || exit 1
echo "Building..."
make || exit 1
if [ -n "$run_tests" ]
then
    echo "Building and running test suite..."
    make check || exit 1
fi
if [ -n "$run_install" ]
then
    echo "Installing..."
    make install || exit 1
fi
cd "$base_dir"
echo "Done!"

