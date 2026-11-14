#include "json.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

namespace climate {
namespace json {

Value Value::from_bool(bool value) {
    Value v;
    v.type = Type::Bool;
    v.boolean = value;
    return v;
}

Value Value::from_number(double value) {
    Value v;
    v.type = Type::Number;
    v.number = value;
    return v;
}

Value Value::from_string(const std::string& value) {
    Value v;
    v.type = Type::String;
    v.text = value;
    return v;
}

Value Value::array() {
    Value v;
    v.type = Type::Array;
    return v;
}

Value Value::object() {
    Value v;
    v.type = Type::Object;
    return v;
}

bool Value::as_bool(bool fallback) const { return type == Type::Bool ? boolean : fallback; }

double Value::as_number(double fallback) const {
    if (type == Type::Number) return number;
    if (type == Type::Bool) return boolean ? 1.0 : 0.0;
    return fallback;
}

std::string Value::as_string(const std::string& fallback) const {
    return type == Type::String ? text : fallback;
}

bool Value::has(const std::string& key) const { return find(key) != nullptr; }

const Value* Value::find(const std::string& key) const {
    if (type != Type::Object) return nullptr;
    for (const auto& member : members)
        if (member.first == key) return &member.second;
    return nullptr;
}

Value& Value::set(const std::string& key, Value value) {
    type = Type::Object;
    for (auto& member : members) {
        if (member.first == key) {
            member.second = value;
            return member.second;
        }
    }
    members.emplace_back(key, std::move(value));
    return members.back().second;
}

Value& Value::push(Value value) {
    type = Type::Array;
    items.push_back(std::move(value));
    return items.back();
}

namespace {

class Parser {
public:
    explicit Parser(const std::string& text) : text_(text) {}

    Value parse_document() {
        skip_whitespace();
        Value value = parse_value();
        skip_whitespace();
        if (position_ != text_.size()) fail("trailing characters after JSON value");
        return value;
    }

private:
    const std::string& text_;
    std::size_t position_ = 0;

    [[noreturn]] void fail(const std::string& message) const {
        std::ostringstream out;
        out << "JSON parse error at offset " << position_ << ": " << message;
        throw std::runtime_error(out.str());
    }

    void skip_whitespace() {
        while (position_ < text_.size()) {
            const char c = text_[position_];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
                ++position_;
            else
                break;
        }
    }

    char peek() const {
        if (position_ >= text_.size()) return '\0';
        return text_[position_];
    }

    char take() {
        if (position_ >= text_.size()) fail("unexpected end of input");
        return text_[position_++];
    }

    void expect(char expected) {
        const char c = take();
        if (c != expected) {
            std::string message = "expected '";
            message.push_back(expected);
            message.push_back('\'');
            fail(message);
        }
    }

    Value parse_value() {
        skip_whitespace();
        const char c = peek();
        if (c == '\0') fail("unexpected end of input");
        if (c == '{') return parse_object();
        if (c == '[') return parse_array();
        if (c == '"') return Value::from_string(parse_string());
        if (c == 't') {
            consume_literal("true");
            return Value::from_bool(true);
        }
        if (c == 'f') {
            consume_literal("false");
            return Value::from_bool(false);
        }
        if (c == 'n') {
            consume_literal("null");
            return Value();
        }
        return parse_number();
    }

    void consume_literal(const char* literal) {
        for (const char* p = literal; *p != '\0'; ++p) {
            if (take() != *p) fail("invalid literal");
        }
    }

    Value parse_object() {
        Value object = Value::object();
        expect('{');
        skip_whitespace();
        if (peek() == '}') {
            ++position_;
            return object;
        }
        while (true) {
            skip_whitespace();
            const std::string key = parse_string();
            skip_whitespace();
            expect(':');
            Value value = parse_value();
            object.members.emplace_back(key, std::move(value));
            skip_whitespace();
            const char c = take();
            if (c == '}') break;
            if (c != ',') fail("expected ',' or '}' in object");
        }
        return object;
    }

    Value parse_array() {
        Value array = Value::array();
        expect('[');
        skip_whitespace();
        if (peek() == ']') {
            ++position_;
            return array;
        }
        while (true) {
            Value value = parse_value();
            array.items.push_back(std::move(value));
            skip_whitespace();
            const char c = take();
            if (c == ']') break;
            if (c != ',') fail("expected ',' or ']' in array");
        }
        return array;
    }

    std::string parse_string() {
        skip_whitespace();
        expect('"');
        std::string result;
        while (true) {
            const char c = take();
            if (c == '"') break;
            if (c != '\\') {
                result.push_back(c);
                continue;
            }
            const char escape_char = take();
            switch (escape_char) {
                case '"': result.push_back('"'); break;
                case '\\': result.push_back('\\'); break;
                case '/': result.push_back('/'); break;
                case 'b': result.push_back('\b'); break;
                case 'f': result.push_back('\f'); break;
                case 'n': result.push_back('\n'); break;
                case 'r': result.push_back('\r'); break;
                case 't': result.push_back('\t'); break;
                case 'u': {
                    unsigned code = 0;
                    for (int i = 0; i < 4; ++i) {
                        const char h = take();
                        code <<= 4;
                        if (h >= '0' && h <= '9')
                            code |= static_cast<unsigned>(h - '0');
                        else if (h >= 'a' && h <= 'f')
                            code |= static_cast<unsigned>(h - 'a' + 10);
                        else if (h >= 'A' && h <= 'F')
                            code |= static_cast<unsigned>(h - 'A' + 10);
                        else
                            fail("invalid unicode escape");
                    }
                    // Encode as UTF-8 (surrogate pairs are merged when present).
                    if (code >= 0xD800 && code <= 0xDBFF) {
                        if (peek() != '\\') fail("unpaired high surrogate");
                        ++position_;
                        expect('u');
                        unsigned low = 0;
                        for (int i = 0; i < 4; ++i) {
                            const char h = take();
                            low <<= 4;
                            if (h >= '0' && h <= '9')
                                low |= static_cast<unsigned>(h - '0');
                            else if (h >= 'a' && h <= 'f')
                                low |= static_cast<unsigned>(h - 'a' + 10);
                            else if (h >= 'A' && h <= 'F')
                                low |= static_cast<unsigned>(h - 'A' + 10);
                            else
                                fail("invalid unicode escape");
                        }
                        code = 0x10000u + ((code - 0xD800u) << 10) + (low - 0xDC00u);
                    }
                    if (code < 0x80) {
                        result.push_back(static_cast<char>(code));
                    } else if (code < 0x800) {
                        result.push_back(static_cast<char>(0xC0 | (code >> 6)));
                        result.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    } else if (code < 0x10000) {
                        result.push_back(static_cast<char>(0xE0 | (code >> 12)));
                        result.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
                        result.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    } else {
                        result.push_back(static_cast<char>(0xF0 | (code >> 18)));
                        result.push_back(static_cast<char>(0x80 | ((code >> 12) & 0x3F)));
                        result.push_back(static_cast<char>(0x80 | ((code >> 6) & 0x3F)));
                        result.push_back(static_cast<char>(0x80 | (code & 0x3F)));
                    }
                    break;
                }
                default: fail("invalid escape sequence");
            }
        }
        return result;
    }

    Value parse_number() {
        const std::size_t start = position_;
        if (peek() == '-' || peek() == '+') ++position_;
        while (position_ < text_.size()) {
            const char c = text_[position_];
            if ((c >= '0' && c <= '9') || c == '.' || c == 'e' || c == 'E' || c == '+' || c == '-')
                ++position_;
            else
                break;
        }
        if (start == position_) fail("expected a value");
        const std::string token = text_.substr(start, position_ - start);
        char* end = nullptr;
        const double value = std::strtod(token.c_str(), &end);
        if (end == nullptr || *end != '\0') fail("invalid number '" + token + "'");
        return Value::from_number(value);
    }
};

void dump_into(const Value& value, std::string& out, int indent, int depth) {
    const auto newline = [&](int level) {
        if (indent < 0) return;
        out.push_back('\n');
        out.append(static_cast<std::size_t>(indent) * static_cast<std::size_t>(level), ' ');
    };
    switch (value.type) {
        case Value::Type::Null: out += "null"; break;
        case Value::Type::Bool: out += value.boolean ? "true" : "false"; break;
        case Value::Type::Number: {
            if (!std::isfinite(value.number)) {
                out += "null";
                break;
            }
            std::ostringstream stream;
            stream.precision(12);
            stream << value.number;
            out += stream.str();
            break;
        }
        case Value::Type::String: out += escape(value.text); break;
        case Value::Type::Array: {
            if (value.items.empty()) {
                out += "[]";
                break;
            }
            out.push_back('[');
            for (std::size_t i = 0; i < value.items.size(); ++i) {
                if (i > 0) out.push_back(',');
                newline(depth + 1);
                dump_into(value.items[i], out, indent, depth + 1);
            }
            newline(depth);
            out.push_back(']');
            break;
        }
        case Value::Type::Object: {
            if (value.members.empty()) {
                out += "{}";
                break;
            }
            out.push_back('{');
            for (std::size_t i = 0; i < value.members.size(); ++i) {
                if (i > 0) out.push_back(',');
                newline(depth + 1);
                out += escape(value.members[i].first);
                out += indent < 0 ? ":" : ": ";
                dump_into(value.members[i].second, out, indent, depth + 1);
            }
            newline(depth);
            out.push_back('}');
            break;
        }
    }
}

}  // namespace

std::string escape(const std::string& text) {
    std::string out;
    out.push_back('"');
    for (char c : text) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\b': out += "\\b"; break;
            case '\f': out += "\\f"; break;
            case '\n': out += "\\n"; break;
            case '\r': out += "\\r"; break;
            case '\t': out += "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buffer[8];
                    std::snprintf(buffer, sizeof(buffer), "\\u%04x",
                                  static_cast<unsigned>(static_cast<unsigned char>(c)));
                    out += buffer;
                } else {
                    out.push_back(c);
                }
        }
    }
    out.push_back('"');
    return out;
}

Value parse(const std::string& text) { return Parser(text).parse_document(); }

std::string dump(const Value& value, int indent) {
    std::string out;
    dump_into(value, out, indent, 0);
    return out;
}

}  // namespace json
}  // namespace climate
