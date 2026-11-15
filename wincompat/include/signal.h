#ifndef CLIMATE_WINCOMPAT_SIGNAL_H
#define CLIMATE_WINCOMPAT_SIGNAL_H
/* Windows has no SIGPIPE; define it so `signal(SIGPIPE, SIG_IGN)` compiles.
 * signal() rejects the unknown number, which is exactly the no-op we want:
 * closed sockets return an error from send() instead of killing the process. */
#include_next <signal.h>
#ifndef SIGPIPE
#define SIGPIPE 13
#endif
#endif
