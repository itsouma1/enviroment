#pragma once

#include <string>
#include <utility>
#include <vector>

namespace climate {
namespace json {

class Value {
public:
    enum class Type { Null, Bool, Number, String, Array, Object };

    Type type = Type::Null;
    bool boolean = false;
    double number = 0.0;
    std::string text;
    std::vector<Value> items;
    std::vector<std::pair<std::string, Value>> members;

    Value() = default;

    static Value from_bool(bool value);
    static Value from_number(double value);
    static Value from_string(const std::string& value);
    static Value array();
    static Value object();

    bool is_null() const { return type == Type::Null; }
    bool is_bool() const { return type == Type::Bool; }
    bool is_number() const { return type == Type::Number; }
    bool is_string() const { return type == Type::String; }
    bool is_array() const { return type == Type::Array; }
    bool is_object() const { return type == Type::Object; }

    bool as_bool(bool fallback) const;
    double as_number(double fallback) const;
    std::string as_string(const std::string& fallback) const;

    bool has(const std::string& key) const;
    const Value* find(const std::string& key) const;
    Value& set(const std::string& key, Value value);
    Value& push(Value value);
};

Value parse(const std::string& text);
std::string dump(const Value& value, int indent = -1);
std::string escape(const std::string& text);

}  // namespace json
}  // namespace climate
