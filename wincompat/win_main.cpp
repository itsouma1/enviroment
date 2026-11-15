/* Windows entry point for the climate engine.
 *
 * The engine's own entry point becomes climate_main at build time
 * (-Dmain=climate_main), so src/ is unchanged. This wrapper exists for two
 * reasons, both about the user seeing what happened:
 *
 *   1. A crash in a console program closes the window before anything can be
 *      read. The unhandled-exception filter prints the exception code and
 *      address, then waits for a key.
 *   2. Double-clicking the exe runs it with no arguments, which prints usage and
 *      exits instantly. In that case the wrapper explains how to start the
 *      workbench and waits, instead of flashing a window.
 */
#include <windows.h>

#include <cstdio>
#include <cstdlib>
#include <exception>

#include <conio.h>

int climate_main(int argc, char** argv);

namespace {

bool launched_without_arguments = false;

void wait_for_key() {
    std::fputs("\nPress any key to close this window...", stdout);
    std::fflush(stdout);
    _getch();
}

LONG WINAPI climate_crash_handler(EXCEPTION_POINTERS* info) {
    const DWORD code = info != nullptr && info->ExceptionRecord != nullptr
                           ? info->ExceptionRecord->ExceptionCode
                           : 0;
    const void* address = info != nullptr && info->ExceptionRecord != nullptr
                              ? info->ExceptionRecord->ExceptionAddress
                              : nullptr;
    std::fflush(stdout);
    std::fprintf(stderr, "\nclimate.exe stopped unexpectedly.\n");
    std::fprintf(stderr, "  exception code: 0x%08lX\n", static_cast<unsigned long>(code));
    std::fprintf(stderr, "  address       : %p\n", address);
    std::fprintf(stderr,
                 "  report this code together with the command you ran; the engine itself did "
                 "not raise a C++ exception.\n");
    std::fflush(stderr);
    wait_for_key();
    return EXCEPTION_EXECUTE_HANDLER;
}

}  // namespace

int main(int argc, char** argv) {
    SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
    SetUnhandledExceptionFilter(climate_crash_handler);
    launched_without_arguments = argc <= 1;

    if (launched_without_arguments) {
        /* Run the engine so it can print its usage text, then stay open. */
        int code = 0;
        try {
            code = climate_main(argc, argv);
        } catch (const std::exception& error) {
            std::fprintf(stderr, "error: %s\n", error.what());
            code = 1;
        } catch (...) {
            std::fprintf(stderr, "error: unknown failure\n");
            code = 1;
        }
        std::fputs(
            "\nThis is a command-line engine. To use the workbench in your browser, either\n"
            "  double-click  \"Start Climate Workbench.bat\"  in this folder, or run\n"
            "  climate.exe ui --port 8080\n"
            "and open http://127.0.0.1:8080 in your browser.\n",
            stdout);
        wait_for_key();
        return code;
    }

    try {
        return climate_main(argc, argv);
    } catch (const std::exception& error) {
        std::fflush(stdout);
        std::fprintf(stderr, "\nclimate.exe failed: %s\n", error.what());
        std::fflush(stderr);
        wait_for_key();
        return 1;
    } catch (...) {
        std::fflush(stdout);
        std::fputs("\nclimate.exe failed with an unknown error.\n", stderr);
        std::fflush(stderr);
        wait_for_key();
        return 1;
    }
}
