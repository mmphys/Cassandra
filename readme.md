# Cassandra

## Pre-requisites

1. C++, Fortran compilers and (optionally) Python, e.g. [Xcode] (or Command Line Tools) and GNU Fortran
2. [ROOT6]
3. [LHAPDF] (I used version 6.2.1) with `NNPDF31_nnlo_as_0118`
4. [APFEL] (I used version 3.0.2 - NB: 3.0.3 didn't work!)

[Xcode]: https://developer.apple.com/xcode/
[ROOT6]: https://root.cern.ch/
[LHAPDF]: https://lhapdf.hepforge.org
[APFEL]: https://apfel.hepforge.org

### Installation on Mac OS

#### 1. [Xcode][Xcode], GNU Fortran and Python

Install [Xcode][Xcode] or the Command Line Tools. 

##### GNU Fortran is part of the GNU Compiler Collection and can be installed from [MacPorts] (I used version 13)

[MacPorts]: https://macports.org

    sudo port install gcc13

Because we will mix `Clang` and `GNU Fortran` we **do not** want to activate `gcc`. We will instead select `GNU Fortran` with `F77=gfortran-mp-13`.  

##### Python (I used version 3.1.2) can also be installed from [MacPorts] and activated

    sudo port install python312
    sudo port select python  python312
    sudo port select python3 python312
    sudo port install py312-cython
    sudo port select cython cython312

Add the following to `PYTHONPATH` (in your `.profile` and current terminal session)

    export PYTHONPATH=${PYTHONPATH:+${PYTHONPATH}:}$HOME/.localMSc/lib/python3.12/site-packages

#### 2. [ROOT6][ROOT6]

I used MacPorts with the cocoa (Mac OS GUI) variant 

    sudo port install root6 +cocoa

#### 3. [LHAPDF][LHAPDF]

Download, then build from source with

    ./configure F77=gfortran-mp-13 --prefix=$HOME/.localMSc
    make
    make install

Install the PDF set used in the thesis `NNPDF31_nnlo_as_0118`

    lhapdf install NNPDF31_nnlo_as_0118

If that failed (e.g. you didn't build LHAPDF with python), you can install this manually

    wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/NNPDF31_nnlo_as_0118.tar.gz -O- | tar xz -C $HOME/.localMSc/share/LHAPDF

#### 4. [APFEL][APFEL]

Download

    git clone https://github.com/scarrazza/apfel.git
    cd apfel

Optionally, select version 3.0.2 with

    git checkout tags/3.0.2

Build from source with

    ./configure F77=gfortran-mp-13 --prefix=$HOME/.localMSc
    make
    make install

Test `APFEL`

    CheckAPFEL

#### 5. Installation troubleshooting

The following shouldn't be necessary on a Mac with build tools already installed, but might be handy to know. 

You can find the `Xcode` SDK directory for Mac OS using

    xcrun --sdk macosx --show-sdk-path 

You can force `clang` and `gcc` to use the latest `Xcode` SDK using

    SDKROOT="$(xcrun --sdk macosx --show-sdk-path)"

To see which SDKs are available (see [stack**overflow**](https://stackoverflow.com/questions/18741675/how-to-get-the-path-of-latest-sdk-available-on-mac) for more)

    xcodebuild -showsdks

[ROOT6] appears to be hardcoded to expect `Command Line Tools` in a specific location. It's a bit clunky, but you can run this (or add it to your .profile)'

    # Root6 is hardcoded to use Command Line Tools here (but shouldn't)
    # If doesn't exist, point at Xcode MacOS platform tools
    (
    CLT=CommandLineTools
    CLTd=/Library/Developer
    if ! [ -e "$CLTd/$CLT" ]; then
      echo "Root6 location for Command Line Tools is hardcoded - fixing"
      set -x
      cd "$CLTd"
      CLTd="$(xcrun --sdk macosx --show-sdk-path)"
      CLTd="${CLTd%/*}"
      CLTd="${CLTd%/*}"
      sudo ln -s "$CLTd" "$CLT"
    fi
    )


# Running a prediction

There are two project executables:

1. **Cassandra** (pass 1) Add F2 predictions to DIS data
2. **Twiggy** (pass 2) Model creation / validation

## **Cassandra** (pass 1)

Run the scripts `B`, `J` and `S` which outputs to `Pass1`

## Create input for pass 2

Run the script `MakePass2.sh` which creates `Pass2` from `Pass1`

### Manually (*deprecated*)

Manually create the Pass2 directory. For each .F2 file created in pass 1, make two copies of the file in Pass2 directory, one for Deuteron data only and one for Proton data only. Update the file names to match the convention for pass 2, i.e.:

* Deuteron_High_BCDMS_F2
* Deuteron_Low_BCDMS_F2
* Deuteron_Low_JLAB_F2
* Proton_High_BCDMS_F2
* Proton_Low_BCDMS_F2
* Proton_Low_JLAB_F2
* NNPDF/Deuteron_High_SLAC_NNPDF_F2
* NNPDF/Deuteron_Low_SLAC_NNPDF_F2
* NNPDF/Proton_High_SLAC_NNPDF_F2
* NNPDF/Proton_Low_SLAC_NNPDF_F2
* reanalyzed/Deuteron_High_SLAC_reanalyzed_F2
* reanalyzed/Deuteron_Low_SLAC_reanalyzed_F2
* reanalyzed/Proton_High_SLAC_reanalyzed_F2
* reanalyzed/Proton_Low_SLAC_reanalyzed_F2

Manually edit the files:

1. Delete the data you don't want (Deuteron or Proton data)
2. Append the target name to the end of the first line in the header, e.g. `Cassandra Prediction BCDMS F2 Proton`
3. Update the Count on line 5 to the number of rows in the file

## **Twiggy** (pass 2)

Edit the script T to choose SLAC_TYPE, DATASET and TARGET, then run the script T

There are plenty of other parameters in the script you can change - once you are up-to-speed, e.g. you can specify a specific model (created in Pass1 at one specific cut) to run against every cut in pass 2 by uncommenting this line:

`#MODEL_NUM="_5"`

NB: There are more options in the scripts, and far more options supported by the executables.
