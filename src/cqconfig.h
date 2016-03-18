#ifndef CQ_CQCONFIG_H
#define CQ_CQCONFIG_H

#define CQ_VERSION_MAJOR 1
#define CQ_VERSION_MINOR 0
#define CQ_VERSION_PATCH 0

#include "timestamp.h"
#include "gitrevision.h"

#define CQ_UTCTIMESTAMP _CMAKE_UTCTIME_
#define CQ_GITREVISION_LONG _CMAKE_GITREVISION_LONG_
#define CQ_GITREVISION_SHORT _CMAKE_GITREVISION_SHORT_

#define CQ_BUILD_YEAR \
    ( \
        (CQ_UTCTIMESTAMP[ 0] - '0') * 1000 +\
        (CQ_UTCTIMESTAMP[ 1] - '0') *  100 +\
        (CQ_UTCTIMESTAMP[ 2] - '0') *   10 +\
        (CQ_UTCTIMESTAMP[ 3] - '0')\
    )

#endif // CQ_CQCONFIG_H

