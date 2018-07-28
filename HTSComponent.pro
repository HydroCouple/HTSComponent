#Author Caleb Amoa Buahin
#Email caleb.buahin@gmail.com
#Date 2018
#License GNU Lesser General Public License (see <http: //www.gnu.org/licenses/> for details).
#The HTSComponent is a stream surface transient storage zone temperature and solute transport model.

TEMPLATE = lib
VERSION = 1.0.0
TARGET = HTSComponent
QT -= gui
QT += testlib

DEFINES += HTSCOMPONENT_LIBRARY
DEFINES += USE_OPENMP
DEFINES += USE_MPI
DEFINES += USE_CVODE
DEFINES += USE_NETCDF
DEFINES += USE_CHPC

#Compile as library or executable
contains(DEFINES,HTSCOMPONENT_LIBRARY){
  TEMPLATE = lib
  message("Compiling HTSComponent as library")
} else {
  TEMPLATE = app
  CONFIG-=app_bundle
  message("Compiling HTSComponent as application")
}

CONFIG += c++11
CONFIG += debug_and_release
CONFIG += optimize_full

PRECOMPILED_HEADER = ./include/stdafx.h

INCLUDEPATH += .\
               ./include \
               ./../HydroCouple/include \
               ./../HydroCoupleSDK/include


HEADERS += ./include/stdafx.h\
           ./include/htscomponent_global.h \
           ./include/htscomponent.h \
           ./include/htscomponentinfo.h \
           ./include/htsmodel.h \
           ./include/iboundarycondition.h \
           ./include/elementjunction.h \
           ./include/odesolver.h \
           ./include/test/htscomponenttest.h \
           ./include/variable.h \
           ./include/element.h \
           ./include/iboundarycondition.h \
           ./include/abstracttimeseriesbc.h \
           ./include/iboundarycondition.h \
           ./include/pointsrctimeseriesbc.h \
           ./include/nonpointsrctimeseriesbc.h \
           ./include/hydraulicstimeseriesbc.h \
           ./include/radiativefluxtimeseriesbc.h \
           ./include/boundarycondition.h \
           ./include/elementinput.h \
           ./include/elementoutput.h


SOURCES +=./src/stdafx.cpp \
          ./src/htscomponent.cpp \
          ./src/htscomponentinfo.cpp \ 
          ./src/main.cpp \
          ./src/odesolver.cpp \
          ./src/element.cpp \
          ./src/htsmodel.cpp \
          ./src/elementjunction.cpp \
          ./src/test/htscomponenttest.cpp \
          ./src/htsmodelio.cpp \
          ./src/htscompute.cpp \
          ./src/abstracttimeseriesbc.cpp \
          ./src/pointsrctimeseriesbc.cpp \
          ./src/nonpointsrctimeseriesbc.cpp \
          ./src/hydraulicstimeseriesbc.cpp \
          ./src/radiativefluxtimeseriesbc.cpp \
          ./src/boundarycondition.cpp \
          ./src/elementinput.cpp \
          ./src/elementoutput.cpp

macx{

    INCLUDEPATH += /usr/local \
                   /usr/local/include

    contains(DEFINES, USE_CVODE){

        message("CVODE enabled")
           LIBS += -L/usr/local/lib -lsundials_cvode
         }

        contains(DEFINES, USE_NETCDF){
        message("NetCDF enabled")
        LIBS += -L/usr/local/lib -lnetcdf-cxx4
    }

    contains(DEFINES,USE_OPENMP){

        QMAKE_CC = /usr/local/opt/llvm/bin/clang
        QMAKE_CXX = /usr/local/opt/llvm/bin/clang++
        QMAKE_LINK = /usr/local/opt/llvm/bin/clang++

        QMAKE_CFLAGS+= -fopenmp
        QMAKE_LFLAGS+= -fopenmp
        QMAKE_CXXFLAGS+= -fopenmp

        INCLUDEPATH += /usr/local/opt/llvm/lib/clang/5.0.0/include
        LIBS += -L /usr/local/opt/llvm/lib -lomp

        message("OpenMP enabled")

     } else {

        message("OpenMP disabled")

     }

    contains(DEFINES,USE_MPI){

        QMAKE_CC = /usr/local/bin/mpicc
        QMAKE_CXX = /usr/local/bin/mpicxx
        QMAKE_LINK = /usr/local/bin/mpicxx

        LIBS += -L/usr/local/lib -lmpi

        message("MPI enabled")

     } else {

        message("MPI disabled")
     }
}

linux{

    INCLUDEPATH += /usr/include \
                   ../gdal/include

    contains(DEFINES,USE_CHPC){

         INCLUDEPATH += /uufs/chpc.utah.edu/sys/installdir/hdf5/1.8.17-c7/include \
                        /uufs/chpc.utah.edu/sys/installdir/netcdf-c/4.3.3.1/include \
                        /uufs/chpc.utah.edu/sys/installdir/netcdf-cxx/4.3.0-c7/include


         LIBS += -L/uufs/chpc.utah.edu/sys/installdir/hdf5/1.8.17-c7/lib -l:libhdf5.so.10.2.0 \
                 -L/uufs/chpc.utah.edu/sys/installdir/netcdf-c/4.4.1/lib -l:libnetcdf.so.11.0.3 \
                 -L/uufs/chpc.utah.edu/sys/installdir/netcdf-cxx/4.3.0-c7/lib -l:libnetcdf_c++4.so.1.0.3

        contains(DEFINES,USE_CVODE){

            message("CVODE enabled")

            INCLUDEPATH += ../sundials-3.1.1/instdir/include
            LIBS += -L../sundials-3.1.1/instdir/lib -lsundials_cvode
        }

         message("Compiling on CHPC")
     }

    contains(DEFINES,USE_OPENMP){

      QMAKE_CFLAGS += -fopenmp
      QMAKE_LFLAGS += -fopenmp
      QMAKE_CXXFLAGS += -fopenmp

      LIBS += -L/usr/lib/x86_64-linux-gnu -lgomp

      message("OpenMP enabled")

      } else {

      message("OpenMP disabled")

     }

    contains(DEFINES,USE_MPI){

        QMAKE_CC = mpicc
        QMAKE_CXX = mpic++
        QMAKE_LINK = mpic++

        QMAKE_CFLAGS += $$system(/usr/local/bin/mpicc --showme:compile)
        QMAKE_CXXFLAGS += $$system(/usr/local/bin/mpic++ --showme:compile)
        QMAKE_LFLAGS += $$system(/usr/local/bin/mpic++ --showme:link)

        LIBS += -L/usr/local/lib/ -lmpi

        message("MPI enabled")

        } else {

        message("MPI disabled")

    }
}

win32{

    #Windows vspkg package manager installation path
    VCPKGDIR = C:/vcpkg/installed/x64-windows

    INCLUDEPATH += $${VCPKGDIR}/include \
                   $${VCPKGDIR}/include/gdal


    CONFIG(debug, debug|release) {
    LIBS += -L$${VCPKGDIR}/debug/lib -lgdald
        } else {
    LIBS += -L$${VCPKGDIR}/lib -lgdal
    }

    contains(DEFINES, USE_CVODE){
    message("CVODE enabled")
    INCLUDEPATH += $${VCPKGDIR}/include/cvode
    CONFIG(release, debug|release) {
    LIBS += -L$${VCPKGDIR}/lib -lsundials_cvode
        } else {
    LIBS += -L$${VCPKGDIR}/debug/lib -lsundials_cvode
        }
    }

    contains(DEFINES, USE_NETCDF){
    message("NetCDF enabled")
    CONFIG(release, debug|release) {
        LIBS += -L$${VCPKGDIR}/lib -lnetcdf \
                -L$${VCPKGDIR}/lib -lnetcdf-cxx4
        } else {
        LIBS += -L$${VCPKGDIR}/debug/lib -lnetcdf \
                -L$${VCPKGDIR}/debug/lib -lnetcdf-cxx4
        }
    }
    
    contains(DEFINES,USE_OPENMP){

        QMAKE_CFLAGS += -openmp
        QMAKE_CXXFLAGS += -openmp

        message("OpenMP enabled")

     } else {

      message("OpenMP disabled")

     }

    contains(DEFINES,USE_MPI){
       message("MPI enabled")

        CONFIG(debug, debug|release) {
            LIBS += -L$${VCPKGDIR}/debug/lib -lmsmpi
          } else {
            LIBS += -L$${VCPKGDIR}/lib -lmsmpi
        }

    } else {
      message("MPI disabled")
    }

    QMAKE_CXXFLAGS += /MP
    QMAKE_LFLAGS += /incremental /debug:fastlink

}

CONFIG(debug, debug|release) {

    win32 {
       QMAKE_CXXFLAGS += /MDd /O2
    }

    macx {
       QMAKE_CXXFLAGS += -O3
    }

    linux {
       QMAKE_CXXFLAGS += -O3
    }

   DESTDIR = ./build/debug
   OBJECTS_DIR = $$DESTDIR/.obj
   MOC_DIR = $$DESTDIR/.moc
   RCC_DIR = $$DESTDIR/.qrc
   UI_DIR = $$DESTDIR/.ui

   macx{
    QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK.* ./build/debug/";
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK
   }

   linux{
    QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/build/debug/*HydroCoupleSDK.* ./build/debug/";
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK
   }

   win32{
    QMAKE_POST_LINK += "copy /B .\..\HydroCoupleSDK\build\debug\HydroCoupleSDK* .\build\debug"
    LIBS += -L./../HydroCoupleSDK/build/debug -lHydroCoupleSDK1
   }
}

CONFIG(release, debug|release) {

   win32 {
     QMAKE_CXXFLAGS += /MD
   }

    RELEASE_EXTRAS = ./build/release
    OBJECTS_DIR = $$RELEASE_EXTRAS/.obj
    MOC_DIR = $$RELEASE_EXTRAS/.moc
    RCC_DIR = $$RELEASE_EXTRAS/.qrc
    UI_DIR = $$RELEASE_EXTRAS/.ui

   macx{
    LIBS += -L./../HydroCoupleSDK/lib/macx -lHydroCoupleSDK
   }

   linux{
    LIBS += -L./../HydroCoupleSDK/lib/linux -lHydroCoupleSDK
   }

   win32{
    LIBS += -L./../HydroCoupleSDK/lib/win32 -lHydroCoupleSDK1
   }

     contains(DEFINES,HTSCOMPONENT_LIBRARY){
         #MacOS
         macx{
             DESTDIR = lib/macx
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/macx/*HydroCoupleSDK.* ./lib/macx/";
          }

         #Linux
         linux{
             DESTDIR = lib/linux
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/linux/*HydroCoupleSDK.* ./lib/linux/";
          }

         #Windows
         win32{
             DESTDIR = lib/win32
             QMAKE_POST_LINK += "copy /B .\..\HydroCoupleSDK\lib\win32\HydroCoupleSDK* .\lib\win32"
          }
     } else {
         #MacOS
         macx{
             DESTDIR = bin/macx
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/macx/*HydroCoupleSDK.* ./bin/macx/";
          }

         #Linux
         linux{
             DESTDIR = bin/linux
             QMAKE_POST_LINK += "cp -a ./../HydroCoupleSDK/lib/linux/*HydroCoupleSDK.* ./bin/linux/";
          }

         #Windows
         win32{
             DESTDIR = bin/win32
             QMAKE_POST_LINK += "copy /B .\..\HydroCoupleSDK\lib\win32\HydroCoupleSDK* .\bin\win32"
          }
     }
}
