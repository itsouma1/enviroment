#ifndef CLIMATE_WINCOMPAT_UNISTD_H
#define CLIMATE_WINCOMPAT_UNISTD_H
#include "win_socket_shim.h"
#include <sys/types.h>
#ifndef _SSIZE_T_DEFINED
#ifndef _SSIZE_T_DEFINED_
typedef long long ssize_t;
#define _SSIZE_T_DEFINED
#endif
#endif
#endif
