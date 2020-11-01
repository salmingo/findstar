bin_PROGRAMS=findstar
findstar_SOURCES=ACatalog.cpp ACatTycho2.cpp ACatUCAC4.cpp AMath.cpp ATimeSpace.cpp findstar.cpp

if DEBUG
  AM_CFLAGS = -g3 -O0 -Wall -DNDEBUG
  AM_CXXFLAGS = -g3 -O0 -Wall -DNDEBUG
else
  AM_CFLAGS = -O3 -Wall
  AM_CXXFLAGS = -O3 -Wall
endif

findstar_LDADD = -lm -lcfitsio -lboost_filesystem-mt  -lboost_date_time-mt
