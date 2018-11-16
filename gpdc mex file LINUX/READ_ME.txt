#--- gpdc

A mex interface to the geopsy (http://geopsy.org/) software, to
replicate the functionality of the gpdc program
(http://www.geopsy.org/wiki/index.php/Gpdc).

The Matlab mex interface replicates much of the geopsy source code found
in the geopsy source file DispersionReader.cpp, but does not implement
all features (e.g. Grid Mode is not implemented). The mex interface also
adds the ability to specify frequency values, rather than just the
minimum and maximum of the range.

Some defaults differ from the gpdc defaults, and these are noted in the
source code.

#-- Usage:

Usage examples:

  out = gpdc(T, Vp, Vs, d)
  out = gpdc(T, Vp, Vs, d, 'nSamples', n, 'minRange', mn, 'maxRange', mx)
  out = gpdc(T, Vp, Vs, d, 'fV', f)

Where required input:

  T  = Thickness(m)    : double vector
  Vp = Vp (m/s)        : double vector
  Vs = Vs (m/s)        : double vector
  d  = density (kg/m3) : double vector

Optional input for frequency values. One or more of these options can be
set:

  n  = number of samples : scalar integer
  mn = minimum range     : scalar double
  mx = maximum range     : scalar double

Alternatively, the frequency values can be set from a vector of doubles:

  f = frequency vector : double vector

If 'fV' is provided 'nSamples', 'minRange' and 'maxRange' values will be
ignored.

Example using fV:

  out = gpdc(T, Vp, Vs, d, 'fV', [10, 20, 30]); 

Which should be equivalent to:

  gpdc -R 5 -s frequency -n 3 -min 10.0 -max 30.0 test.model

#- Output:

  out = 2d matrix

    col1 : x values
    col2 : y values for mode 0
    col3 : y values for mode 1
    ...

For example:

  T  = [7.5, 25, 0];
  Vp = [500, 1350, 2000];
  Vs = [200, 210, 1000];
  d  = [1700, 1900, 2500];
  out = gpdc(T, Vp, Vs, d);
  plot(out(:, 1), out(:, 2:end))

Which should be equivalent to (these command line options match our defaults):

  gpdc -R 5 -s frequency -n 150 -min 1.0 -max 150.0 test.model | figue -c

#-- Compiling:

The mex interface has been mostly tested on CentOS 7, and Windows 7,
with Matlab 2017a and geopsy 2.10.0.

#- Linux:

Presuming the Matlab bin directory is in the PATH, the geopsy software
is installed at ${GEOPSY_HOME}, and we are using the system version of
Qt, installed within the /usr directory, the mex file can be compiled
with:

  mex gpdc.cpp \
    -I/usr/include/QtCore \
    -L${GEOPSY_HOME} \
    -lQGpCoreTools \
    -lQGpCoreWave \
    LDFLAGS="-Wl,-rpath,${GEOPSY_HOME}/lib"

Which should produce the file:

  gpdc.mexa64

Testing has found that if the system version of Qt differs from the
version of Qt used by Matlab, this can cause crashes when using the
Matlab graphical interface.

To avoid the crashes, the geopsy code and mex code can be compiled using
the same version of Qt as used by Matlab. For example, Matlab 2017a uses
Qt 5.5.1. The build process is a bit more complicated, but does not
require building the entire geopsy suite, as only QGpCoreTools and
QGpCoreWave are used.

Notes on building Qt 5.5.1, required geopsy libraries, and mex interface
for Matlab 2017a, where everything will be built and installed within
the current working directory. Presumes Matlab is installed at
${MATLAB_HOME}:

  export CFLAGS='-O2 -fPIC'
  export CXXFLAGS='-O2 -fPIC'
  export FFLAGS='-O2 -fPIC'
  export FCFLAGS='-O2 -fPIC'

  unset QTDIR
  unset QTINC
  unset QTLIB

  #- Qt 5.5.1:

  wget "http://download.qt.io/archive/qt/5.5/5.5.1/single/qt-everywhere-opensource-src-5.5.1.tar.xz"
  tar xJf qt-everywhere-opensource-src-5.5.1.tar.xz

  mkdir build.qt-everywhere-opensource-src-5.5.1
  cd build.qt-everywhere-opensource-src-5.5.1

  mkdir ../qt5
  qt5Prefix="$(readlink -f ../qt5)"

  ../qt-everywhere-opensource-src-5.5.1/configure \
    -opensource \
    -confirm-license \
    -no-rpath \
    -qt-xcb \
    -no-audio-backend \
    -nomake examples \
    -nomake tests \
    -make libs \
    -prefix ${qt5Prefix} \
    -bindir ${qt5Prefix}/bin && \
    LD_LIBRARY_PATH="$(pwd)/lib:${LD_LIBRARY_PATH}" \
    make -j32 && \
    LD_LIBRARY_PATH="$(pwd)/lib:${LD_LIBRARY_PATH}" \
    make -j32 install && \
    cd ..

  PATH="$(pwd)/qt5/bin:${PATH}"
  CPATH="$(pwd)/qt5/include:${CPATH}"
  LIBRARY_PATH="$(pwd)/qt5/lib:${LIBRARY_PATH}"
  LD_LIBRARY_PATH="$(pwd)/qt5/lib:${LIBRARY_PATH}"

  #- geopsy:
  #  sources downloaded from:
  #    http://www.geopsy.org/download.php

  tar xzf geopsypack-55items-src-2.10.1.tar.gz
  cd geopsypack-55items-src-2.10.1

  mkdir ../geopsy
  gpPrefix="$(readlink -f ../geopsy)"

  echo 'yes' | \
    ./configure \
    -prefix ${gpPrefix} \
    -I ${qt5Prefix}/include \
    -I ${qt5Prefix}/include/QtCore \
    -L ${qt5Prefix}/lib \
    -I ${MATLAB_HOME}/extern/include \
    -L ${MATLAB_HOME}/bin/glnxa64

  # Only need to make QGpCoreTools and QGpCoreWave:

  cd QGpCoreTools

  # Patches for Qt 5.5.1:

  cp src/AbstractStream.h \
    src/AbstractStream.h.original
  sed -i 's|fromAscii|fromLatin1|g' \
    src/AbstractStream.h

  cp src/Cache.cpp \
    src/Cache.cpp.original
  sed -i 's|toAscii|toLatin1|g' \
    src/Cache.cpp

  cp src/ConsoleProgress.cpp \
    src/ConsoleProgress.cpp.original
  sed -i 's|toAscii|toLatin1|g' \
    src/ConsoleProgress.cpp

  cp src/CoreApplicationPrivate.cpp \
    src/CoreApplicationPrivate.cpp.original
  sed -i 's|toAscii|toLatin1|g' \
    src/CoreApplicationPrivate.cpp
  sed -i '/qInstallMsgHandler/d' \
    src/CoreApplicationPrivate.cpp

  cp src/XMLClass.cpp \
    src/XMLClass.cpp.original
  sed -i 's|toAscii|toLatin1|g' \
    src/XMLClass.cpp

  cp src/Tar.cpp \
    src/Tar.cpp.original
  sed -i 's|toAscii|toLatin1|g' \
    src/Tar.cpp

  cp src/XMLParser.cpp \
    src/XMLParser.cpp.original
  sed -i 's|toAscii|toLatin1|g' \
    src/XMLParser.cpp

  cp src/Thread.cpp \
    src/Thread.cpp.original
  sed -i '/QThread::finished()/d' \
    src/Thread.cpp

  make -j16 release
  make -j16 release-install
  cd ..

  \cp include/QGpCoreTools/QGpCoreToolsInstallPath.h \
    ${gpPrefix}/include/

  cd QGpCoreWave

  make -j16 release
  make -j16 release-install
  cd ..

  \cp include/QGpCoreWave/QGpCoreWaveInstallPath.h \
    ${gpPrefix}/include/

  cd ..

  #- mex file - presume mex source is already in 'mex' directory ... :

  cd mex

  g++ \
    gpdc.cpp \
    -O2 \
    -fPIC \
    -fpermissive \
    -shared \
    -DMATLAB_MEX_FILE \
    -o gpdc.mexa64 \
    -I${gpPrefix}/include \
    -I${qt5Prefix}/include \
    -I${qt5Prefix}/include/QtCore \
    -I${MATLAB_HOME}/extern/include \
    -L${MATLAB_HOME}/bin/glnxa64 \
    -lmex \
    -lmx \
    -leng \
    -lmat \
    -L${qt5Prefix}/lib \
    -lQt5Core \
    -L${gpPrefix}/lib \
    -lQGpCoreTools \
    -lQGpCoreWave

This will build the file gpdc.mexa64, which will require the libraries
${gpPrefix}/libQGpCoreTools.so.1 and ${gpPrefix}/libQGpCoreWave.so.1 to
run.

For the end user, using the patchelf (https://nixos.org/patchelf.html)
program to set the appropriate RPATH for these files is recommended.

For example, if the mex file was to be installed at and run from the
directory ${HOME}/matlab/gpc, the required files can be copied in to
place:

  ${HOME}/matlab/gpdc/
  ├── gpdc.mexa64
  └── lib/
      ├── libQGpCoreTools.so.1
      └── libQGpCoreWave.so.1

The following patchelf commands can then be used to set appropriate
RPATH values for the files:

  patchelf --remove-rpath ${HOME}/matlab/gpdc/gpdc.mexa64
  patchelf --remove-rpath ${HOME}/matlab/gpdc/lib/libQGpCoreTools.so.1
  patchelf --remove-rpath ${HOME}/matlab/gpdc/lib/libQGpCoreWave.so.1
  patchelf \
    --set-rpath \
    ${MATLAB_HOME}/sys/os/glnxa64:${MATLAB_HOME}/bin/glnxa64:${HOME}/matlab/gpdc/lib \
    gpdc.mexa64

The directory ${HOME}/matlab/gpdc can then be added to the Matlab path
using the preferred method.

#- Windows:

People who are more used to doing this kind of thing in Windows may know
better ways of doing this ...

Building on Windows (Windows 7 / Matlab 2017a) was done within a Git
Bash shell (https://git-scm.com/download/win).

TDM mingw 4.9.2 64 bit compilers were used, downloaded from:

  https://sourceforge.net/projects/tdm-gcc/files/TDM-GCC%20Installer/Previous/1.1309.0/

This installs to the directory:

  C:\TDM-GCC-64

The Qt 4.8.5 / mingw 4.8.2 build was downloaded from:

  https://sourceforge.net/projects/mingwbuilds/files/external-binary-packages/Qt-Builds/

The files were extracted, and moved in to the directory:

  C:\Test

Which seems to be the expected location.

The mex file was then built within the Git Bash shell:

  #- Set variables:

  export QTBASEDIR="/c/Test"
  export QTDIR="${QTBASEDIR}/Qt-4.8.5-x86_64"
  export PATH="${QTDIR}/bin:${PATH}"
  export CPATH="${QTDIR}/include/QtCore:${CPATH}"
  export CPATH="${QTBASEDIR}/prerequisites-x86_64/include:${QTDIR}/include:${CPATH}"
  export LIBRARY_PATH="${QTDIR}/lib:${QTBASEDIR}/mingw64/x86_64-w64-mingw32/lib:${LIBRARY_PATH}"
  export QMAKESPEC="${QTDIR}/mkspecs/win32-g++"
  export MINGWROOT="/c/TDM-GCC-64"
  export PATH="${MINGWROOT}/bin:${PATH}"
  export GEOPSY_HOME="/c/Users/user/geopsy"
  export CPATH="${GEOPSY_HOME}/include:${CPATH}"
  export LIBRARY_PATH="${GEOPSY_HOME}/lib:${LIBRARY_PATH}"
  export CFLAGS='-O2 -fPIC'
  export CXXFLAGS='-O2 -fPIC'
  export FFLAGS='-O2 -fPIC'
  export FCFLAGS='-O2 -fPIC'

  #- Edit ${QMAKESPEC}/qmake.conf:

  cp ${QMAKESPEC}/qmake.conf \
    ${QMAKESPEC}/qmake.conf.original
  sed -i 's|^\(QMAKE_CFLAGS_WARN_ON\).*$|\1    =|g' \
    ${QMAKESPEC}/qmake.conf
  sed -i 's|^\(QMAKE_CFLAGS_RELEASE.*$\)|\1 -fpermissive|g' \
    ${QMAKESPEC}/qmake.conf

  #- Copy ${QTDIR}/lib/QTCore4.dll.bak ${QTDIR}/lib/QTCore.dll

  cp ${QTDIR}/lib/QTCore4.dll.bak \
    ${QTDIR}/lib/QTCore.dll

  #- Extract geopsy source, and configure:

  tar xzf geopsypack-55items-src-2.10.1.tar.gz
  cd geopsy-2.10.1

  ./configure \
    -prefix ${GEOPSY_HOME}

  #  only need to make QGpCoreTools and QGpCoreWave:

  cd QGpCoreTools
  cp QGpCoreTools.pro QGpCoreTools.pro.original
  sed -i 's|^\(DEFINES.*$\)|\1 QT_DISABLE_DEPRECATED_BEFORE=0 Q_WS_WIN=1|g' \
    QGpCoreTools.pro
  sed -i 's|\(LIBS\s.*$\)|\1 -lz|g' \
    QGpCoreTools.pro

  mingw32-make -j16 release
  mingw32-make -j16 release-install
  cd ..

  \cp include/QGpCoreTools/QGpCoreToolsInstallPath.h \
    ${GEOPSY_HOME}/include/

  cd QGpCoreWave
  cp QGpCoreWave.pro QGpCoreWave.pro.original
  sed -i 's|^\(DEFINES.*$\)|\1 QT_DISABLE_DEPRECATED_BEFORE=0 Q_WS_WIN=1|g' \
    QGpCoreWave.pro
  sed -i 's|\(LIBS\s.*$\)|\1 -lz|g' \
    QGpCoreWave.pro

  mingw32-make -j16 release
  mingw32-make -j16 release-install
  cd ..

  \cp include/QGpCoreWave/QGpCoreWaveInstallPath.h \
    ${GEOPSY_HOME}/include/

  cd ..

  #- Compile mex file:

  MATLAB_HOME='/c/Program Files/MATLAB/R2017a'

  g++ \
    gpdc.cpp \
    -O2 \
    -fpermissive \
    -shared \
    -DMATLAB_MEX_FILE \
    -o gpdc.mexw64 \
    -I"${MATLAB_HOME}/extern/include" \
    "${MATLAB_HOME}/extern/lib/win64/mingw64/libmex.lib" \
    "${MATLAB_HOME}/extern/lib/win64/mingw64/libmx.lib" \
    "${MATLAB_HOME}/extern/lib/win64/mingw64/libeng.lib" \
    "${MATLAB_HOME}/extern/lib/win64/mingw64/libmat.lib" \
    -lQGpCoreTools1 \
    -lQGpCoreWave1 \
    -lQtCore4

  #- As well as the mex file, the following libraries are required to
  #  be in the same directory:

  \cp \
    ${GEOPSY_HOME}/lib/*.dll \
    ${QTDIR}/lib/QtCore4.dll.bak \
    ${QTBASEDIR}/mingw64/bin/libwinpthread-1.dll \
    ${QTBASEDIR}/mingw64/bin/libgcc_s_seh-1.dll \
    ${QTBASEDIR}/mingw64/bin/libstdc++-6.dll \
    .

  \mv QtCore4.dll.bak QtCore4.dll

#-

