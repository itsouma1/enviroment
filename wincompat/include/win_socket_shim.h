/* Windows compatibility layer for the POSIX socket code in src/ui.cpp.
 *
 * Nothing in src/ is modified: this header is pulled in through the shim
 * headers arpa/inet.h, netinet/in.h, sys/socket.h and unistd.h, which the
 * build puts ahead of the system include path with -I.
 *
 * Winsock SOCKET handles are 64-bit on x64, while the original code stores
 * them in `int`. To keep that code untouched and still correct, every handle
 * is registered in a small table and the POSIX-style call returns a 32-bit
 * descriptor index that is translated back on each call.
 *
 * `close` is a real function rather than a macro: standard headers in the
 * same translation unit call std::basic_filebuf::close(), and a macro would
 * rewrite that member call.
 */
#ifndef CLIMATE_WIN_SOCKET_SHIM_H
#define CLIMATE_WIN_SOCKET_SHIM_H

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <winsock2.h>
#include <ws2tcpip.h>
#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

int climate_win_socket(int domain, int type, int protocol);
int climate_win_bind(int descriptor, const struct sockaddr* address, int length);
int climate_win_listen(int descriptor, int backlog);
int climate_win_accept(int descriptor, struct sockaddr* address, socklen_t* length);
int climate_win_recv(int descriptor, void* buffer, unsigned int length, int flags);
int climate_win_send(int descriptor, const void* buffer, unsigned int length, int flags);
int climate_win_setsockopt(int descriptor, int level, int option, const void* value, int length);
int close(int descriptor);

#ifdef __cplusplus
}
#endif

#define socket(domain, type, protocol) climate_win_socket(domain, type, protocol)
#define bind(descriptor, address, length) climate_win_bind(descriptor, address, length)
#define listen(descriptor, backlog) climate_win_listen(descriptor, backlog)
#define accept(descriptor, address, length) climate_win_accept(descriptor, address, length)
#define recv(descriptor, buffer, length, flags) climate_win_recv(descriptor, buffer, length, flags)
#define send(descriptor, buffer, length, flags) climate_win_send(descriptor, buffer, length, flags)
#define setsockopt(descriptor, level, option, value, length) \
    climate_win_setsockopt(descriptor, level, option, value, length)

#endif /* CLIMATE_WIN_SOCKET_SHIM_H */
