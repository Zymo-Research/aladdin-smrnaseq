def setup_channel(path, name, required, missing_message) {
    if (path) {
        ch = Channel
                 .fromPath(path, checkIfExists: true)
                 .ifEmpty { exit 1, "$name not found at: $path" }
    } else if (required) {
        exit 1, "$name is a required input!"
    } else {
        ch = Channel.empty()
        log.info "No $name provided - $missing_message"
    }
    return ch
}