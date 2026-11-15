#pragma once

#include <string>

namespace climate {

struct UiOptions {
    std::string host = "127.0.0.1";
    int port = 8080;
};

// Serves the embedded single-page workbench and its JSON API. Blocks until the
// process is interrupted; one short-lived thread handles each request.
int run_ui_server(const UiOptions& options);

}  // namespace climate
