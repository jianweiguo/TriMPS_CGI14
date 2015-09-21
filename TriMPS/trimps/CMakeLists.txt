
IF(WIN32)
	# Od : no optimisation; Ot: favor code speed over code size.
	if(DO_NOT_OPTIMIZE_WINDOWS_RELEASE)
		SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Od /Ot")
	endif(DO_NOT_OPTIMIZE_WINDOWS_RELEASE)
ELSE(WIN32)
	# CGAL needs rounding math to compile well with GCC
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
ENDIF(WIN32)

AUX_SOURCE_DIRECTORY(. SUB_SOURCES)
AUX_SOURCE_DIRECTORY(generated SOURCES)
SET(SOURCES ${SOURCES} ${SUB_SOURCES})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/trimps)
file(GLOB HeaderFiles *.h)

SET(USE_OPENCV false)
IF (USE_OPENCV) 
   FIND_PACKAGE(OPENCV)
   IF(OPENCV_FOUND)
       MESSAGE("OpenCV found ${OPENCV_FILE_NAME}")
   ELSE(OPENCV_FOUND)
	   MESSAGE("OpenCV not found")
   ENDIF(OPENCV_FOUND)
ENDIF (USE_OPENCV) 

IF (WIN32)
   INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/glut-win32/include)
   LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/glut-win32/lib)
   ADD_EXECUTABLE(trimps ${SOURCES} ${HeaderFiles})
ELSE (WIN32)
   LINK_DIRECTORIES(${CGAL_DIR}/lib/)
   LINK_DIRECTORIES(/usr/X11R6/lib)
   ADD_EXECUTABLE(trimps ${SOURCES} ${HeaderFiles})
ENDIF (WIN32)

ADD_DEPENDENCIES(trimps geex)

ADD_DEFINITIONS(-DNOMINMAX) 

IF (WIN32)
   TARGET_LINK_LIBRARIES(trimps geex opengl32 glu32)
   IF(USE_OPENCV)
		TARGET_LINK_LIBRARIES(trimps debug opencv_core231d debug opencv_imgproc231d debug opencv_highgui231d debug opencv_video231d)
		TARGET_LINK_LIBRARIES(trimps optimized opencv_core231 optimized opencv_imgproc231 optimized opencv_highgui231 optimized opencv_video231)
   ENDIF(USE_OPENCV) 
   # SET_TARGET_PROPERTIES(gx_sms2d PROPERTIES LINK_FLAGS "/FORCE:MULTIPLE")
ELSE (WIN32)
   TARGET_LINK_LIBRARIES(trimps geex glut GL GLU CGAL)
ENDIF (WIN32)

