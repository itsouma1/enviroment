/* Implementation of the Windows socket compatibility layer. See the header. */
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <winsock2.h>
#include <ws2tcpip.h>

#include <stddef.h>

#define CLIMATE_WIN_MAX_SOCKETS 256

static SOCKET climate_win_handle[CLIMATE_WIN_MAX_SOCKETS];
static int climate_win_in_use[CLIMATE_WIN_MAX_SOCKETS];
static int climate_win_started = 0;

static void climate_win_startup(void) {
    if (!climate_win_started) {
        WSADATA data;
        WSAStartup(MAKEWORD(2, 2), &data);
        climate_win_started = 1;
    }
}

static int climate_win_register(SOCKET handle) {
    int i;
    if (handle == INVALID_SOCKET) return -1;
    climate_win_startup();
    for (i = 0; i < CLIMATE_WIN_MAX_SOCKETS; ++i) {
        if (!climate_win_in_use[i]) {
            climate_win_in_use[i] = 1;
            climate_win_handle[i] = handle;
            return i;
        }
    }
    closesocket(handle);
    return -1;
}

static SOCKET climate_win_lookup(int descriptor) {
    if (descriptor < 0 || descriptor >= CLIMATE_WIN_MAX_SOCKETS) return INVALID_SOCKET;
    if (!climate_win_in_use[descriptor]) return INVALID_SOCKET;
    return climate_win_handle[descriptor];
}

extern "C" int climate_win_socket(int domain, int type, int protocol) {
    climate_win_startup();
    return climate_win_register(socket(domain, type, protocol));
}

extern "C" int climate_win_bind(int descriptor, const struct sockaddr* address, int length) {
    return bind(climate_win_lookup(descriptor), address, length);
}

extern "C" int climate_win_listen(int descriptor, int backlog) {
    return listen(climate_win_lookup(descriptor), backlog);
}

extern "C" int climate_win_accept(int descriptor, struct sockaddr* address, socklen_t* length) {
    int native_length = length != NULL ? *length : 0;
    const int accepted = accept(climate_win_lookup(descriptor), address, &native_length);
    if (length != NULL) *length = (socklen_t)native_length;
    return climate_win_register(accepted);
}

extern "C" int climate_win_recv(int descriptor, void* buffer, unsigned int length, int flags) {
    return recv(climate_win_lookup(descriptor), (char*)buffer, (int)length, flags);
}

extern "C" int climate_win_send(int descriptor, const void* buffer, unsigned int length, int flags) {
    return send(climate_win_lookup(descriptor), (const char*)buffer, (int)length, flags);
}

extern "C" int climate_win_setsockopt(int descriptor, int level, int option, const void* value, int length) {
    return setsockopt(climate_win_lookup(descriptor), level, option, (const char*)value, length);
}

extern "C" int close(int descriptor) {
    const SOCKET handle = climate_win_lookup(descriptor);
    if (handle == INVALID_SOCKET) return -1;
    climate_win_in_use[descriptor] = 0;
    return closesocket(handle);
}
