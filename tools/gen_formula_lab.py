#!/usr/bin/env python3
"""Generate src/formula_lab_gen.cpp from the public declarations in src/*.hpp.

Every function declared in a public header becomes a runnable formula-lab entry:
the generator emits a call site, the conversions from the shared lab inputs to
the declared parameter types, and a formatter for the return value. Nothing is
hand-picked; re-running the generator after adding a declaration adds a lab
entry.

Run: python3 tools/gen_formula_lab.py
"""

import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def strip_comments(text):
    out = []
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        nxt = text[i + 1] if i + 1 < n else ''
        if c == '/' and nxt == '/':
            while i < n and text[i] != '\n':
                i += 1
            continue
        if c == '/' and nxt == '*':
            i += 2
            while i < n - 1 and not (text[i] == '*' and text[i + 1] == '/'):
                i += 1
            i += 2
            continue
        out.append(c)
        i += 1
    return ''.join(out)


DECL = re.compile(
    r'([A-Za-z_][\w:<>,\.\s\*&]*?)\s*\b([A-Za-z_]\w*)\s*\(([^;{}]*)\)\s*(const)?\s*([;{=])')
STRUCT = re.compile(r'\bstruct\s+([A-Za-z_]\w*)\s*\{')
ENUM = re.compile(r'\benum\s+(?:class\s+)?([A-Za-z_]\w*)\s*(?::\s*[\w:]+)?\s*\{([^{}]*)\}')
ALIAS = re.compile(r'\busing\s+([A-Za-z_]\w*)\s*=\s*([^;]+);')
TYPEDEF = re.compile(r'\btypedef\s+([^;]+?)\s+([A-Za-z_]\w*)\s*;')

SKIP_PREFIX = ('namespace', 'using', 'typedef', 'template', 'return', 'if', 'for',
               'while', 'switch', 'static_assert', 'friend', 'class', 'enum',
               'static', 'inline', 'extern')

STAGES = [
    ('Observation and retrieval', [
        'aerosol', 'cryo', 'earth', 'fetch', 'land', 'netcdf', 'ocean', 'ocean2',
        'ozone', 'precip', 'readers', 'retrieval', 'sstunc']),
    ('Time series and change', [
        'analysis', 'change', 'detect', 'ensemble', 'harmonic', 'hmm', 'infotheory',
        'intercal', 'loess', 'magnitude', 'maturity', 'multitrend', 'quantile',
        'regime', 'series', 'ssa', 'stats', 'stl', 'trend']),
    ('Space, attribution and dynamics', [
        'atmos', 'attrib', 'attrib2', 'causal', 'causal2', 'clustering', 'coherence',
        'copula', 'flow', 'geostat', 'lagrangian', 'landuse', 'pde', 'regional',
        'relevance', 'restrand', 'spatial', 'tempac']),
    ('Uncertainty and validation', [
        'assim', 'climateidx', 'homogen', 'metrology', 'power', 'uncertainty',
        'validation']),
    ('Physics and radiative transfer', ['physics']),
    ('Data, numerics and reasoning', [
        'catalog', 'csv', 'fft', 'intent', 'json', 'kernel', 'linalg', 'ml',
        'report', 'ui']),
]


def stage_of(module):
    name = module[:-4] if module.endswith('.hpp') else module
    for index, (label, modules) in enumerate(STAGES, start=1):
        if name in modules:
            return '%02d' % index, label
    return '07', 'Supporting modules'


def split_top(text):
    parts = []
    depth = 0
    current = ''
    for c in text:
        if c in '(<[':
            depth += 1
        if c in ')>]':
            depth -= 1
        if c == ',' and depth == 0:
            parts.append(current.strip())
            current = ''
        else:
            current += c
    parts.append(current.strip())
    return [p for p in parts if p]


PARAM = re.compile(r'^(?P<type>.*?)(?P<name>[A-Za-z_]\w*)\s*(?:=\s*(?P<default>.+))?$', re.S)


def parse_param(text):
    text = ' '.join(text.split())
    default = None
    if '=' in text:
        head, _, tail = text.partition('=')
        text = head.strip()
        default = tail.strip()
    m = PARAM.match(text)
    if not m:
        return None
    ptype = ' '.join(m.group('type').split())
    if not ptype:
        return None
    return {'type': ptype, 'name': m.group('name'), 'default': default, 'raw': text}


def normalise(ptype):
    t = ptype.replace('const ', '').replace('&', '').replace('*', '').strip()
    return ' '.join(t.split())


SCALARS = {
    'double': 'double', 'float': 'float', 'long double': 'long double',
    'int': 'int', 'long': 'long', 'long long': 'long long', 'short': 'short',
    'unsigned': 'unsigned', 'unsigned int': 'unsigned', 'unsigned long': 'unsigned long',
    'std::size_t': 'std::size_t', 'size_t': 'std::size_t',
    'std::uint64_t': 'std::uint64_t', 'std::uint32_t': 'std::uint32_t',
    'std::int64_t': 'std::int64_t', 'bool': 'bool', 'char': 'char',
}

INT_TYPES = ('int', 'long', 'long long', 'short', 'unsigned', 'unsigned int',
             'unsigned long', 'std::size_t', 'size_t', 'std::uint64_t',
             'std::uint32_t', 'std::int64_t', 'char')

DEFAULTS = [
    ('refractive_imaginary', 0.0), ('imaginary', 0.0),
    ('refractive_real', 1.5), ('poisson', 0.3), ('damage', 0.5), ('fraction', 0.5), ('albedo', 0.3),
    ('concentration', 0.6), ('coverage', 0.5), ('correlation', 0.5),
    ('ratio', 0.5), ('efficiency', 0.5), ('reflectance', 0.2),
    ('importance', 0.5), ('moisture', 0.35), ('buttressing', 0.7),
    ('stress', 0.5), ('strain', 0.05), ('curvature', 0.01),
    ('latitude', 45.0), ('lon', 0.0), ('lat', 45.0),
    ('alpha', 0.05), ('threshold', 0.0), ('level', 0.05),
    ('ppm', 420.0), ('co2', 420.0), ('reference', 280.0),
    ('temperature', 288.15), ('temp', 288.15), ('bt11', 300.2), ('bt12', 299.3),
    ('pressure', 101325.0), ('density', 1025.0), ('rho', 1025.0),
    ('radius', 6371000.0), ('wavelength', 5.5e-7), ('period', 12.0),
    ('window', 12.0), ('lag', 1.0), ('order', 2.0), ('degree', 2.0),
    ('harmonic', 4.0), ('iteration', 100.0), ('seed', 42.0),
    ('confidence', 0.95), ('percentile', 0.5), ('probability', 0.95),
    ('rate', 1.0), ('depth', 10.0), ('height', 1000.0), ('mass', 1.0),
    ('flux', 3.7), ('forcing', 3.7), ('feedback', -1.25),
    ('size', 12.0), ('count', 12.0), ('n', 12.0), ('m', 4.0), ('rows', 4.0),
    ('steps', 32.0), ('sigma', 0.5), ('asymmetry', 0.3), ('step_size', 1.0),
    ('columns', 3.0), ('cols', 3.0), ('years', 24.0), ('months', 120.0),
    ('x', 1.0), ('y', 1.0), ('z', 0.0), ('u10', 10.0), ('v10', 0.0),
    ('salinity', 35.0), ('gravity', 9.80665), ('time_index', 0.0),
    ('scale_index', 1.0), ('index', 0.0), ('steps', 32.0), ('points', 16.0),
    ('tolerance', 1e-6), ('min', 0.0), ('max', 1.0),
]


def default_for(param):
    if param['default'] is not None:
        raw = param['default']
        if re.match(r'^-?\d+(\.\d*)?([eE][-+]?\d+)?$', raw):
            return float(raw)
        if raw in ('true', 'false'):
            return 1.0 if raw == 'true' else 0.0
        return 1.0
    name = param['name'].lower()
    tokens = set(name.split('_'))
    for hint, value in sorted(DEFAULTS, key=lambda item: -len(item[0])):
        if hint == name or hint in tokens:
            return value
    return 1.0


def brace_body(source, open_index):
    """Return the text inside the braces that open at open_index."""
    depth = 0
    for i in range(open_index, len(source)):
        if source[i] == '{':
            depth += 1
        elif source[i] == '}':
            depth -= 1
            if depth == 0:
                return source[open_index + 1:i]
    return source[open_index + 1:]


def scan():
    functions = []
    structs = {}
    enums = {}
    aliases = {}
    modules = []
    seen = set()
    for path in sorted(glob.glob(os.path.join(ROOT, 'src', '*.hpp'))):
        module = os.path.basename(path)
        if module in ('formula_lab.hpp', 'ui_page.hpp'):
            continue
        source = strip_comments(open(path, encoding='utf-8').read())
        modules.append(module)
        for m in ALIAS.finditer(source):
            aliases.setdefault(m.group(1), ' '.join(m.group(2).split()))
        for m in TYPEDEF.finditer(source):
            aliases.setdefault(m.group(2), ' '.join(m.group(1).split()))
        for m in STRUCT.finditer(source):
            members = []
            body = brace_body(source, m.end() - 1)
            for statement in body.split(';'):
                statement = ' '.join(statement.split())
                if not statement or '(' in statement or statement.startswith('using'):
                    continue
                pm = PARAM.match(statement)
                if pm and pm.group('type').strip():
                    members.append({'type': ' '.join(pm.group('type').split()),
                                    'name': pm.group('name')})
            structs.setdefault(m.group(1), members)
        for m in ENUM.finditer(source):
            names = []
            for item in m.group(2).split(','):
                item = item.split('=')[0].strip()
                if re.match(r'^[A-Za-z_]\w*$', item):
                    names.append(item)
            enums.setdefault(m.group(1), names)
        masked = list(source)
        for m in re.finditer(r'\b(?:struct|class|enum)\s+[A-Za-z_]\w*[^;{]*\{', source):
            body = brace_body(source, m.end() - 1)
            start = m.end()
            for offset in range(start, start + len(body) + 1):
                if offset < len(masked):
                    masked[offset] = ' '
        free_source = ''.join(masked)
        for m in DECL.finditer(free_source):
            ret, name, params, _const, term = m.groups()
            if term in ('{', '='):
                continue
            ret = ' '.join(ret.split())
            if not ret or ret.split()[0] in SKIP_PREFIX:
                continue
            parsed = [parse_param(p) for p in split_top(params)]
            if any(p is None for p in parsed):
                continue
            key = (name, tuple(p['raw'] for p in parsed))
            if key in seen:
                continue
            seen.add(key)
            functions.append({'module': module, 'ret': ret, 'name': name,
                              'params': parsed})
    return functions, structs, enums, aliases, modules


class Unsupported(Exception):
    pass


def resolve(ptype, aliases, depth=0):
    """Resolve aliases one level at a time (Matrix, Weights, ComplexField...)."""
    t = normalise(ptype)
    if depth > 6:
        return t
    if t in aliases:
        return resolve(aliases[t], aliases, depth + 1)
    inner = re.match(r'^(std::vector|std::function)\s*<(.*)>$', t)
    if inner:
        kind = inner.group(1)
        items = [resolve(x, aliases, depth + 1) for x in split_top(inner.group(2))]
        return '%s<%s>' % (kind, ', '.join(items))
    pair = re.match(r'^(std::pair)\s*<(.*)>$', t)
    if pair:
        items = [resolve(x, aliases, depth + 1) for x in split_top(pair.group(2))]
        return 'std::pair<%s>' % ', '.join(items)
    return t


def kind_of(ptype, structs, enums, aliases):
    t = resolve(ptype, aliases)
    if t in SCALARS:
        return t, 'number'
    if t in ('std::complex<double>', 'std::complex<float>'):
        return t, 'complex'
    if t in ('std::vector<double>', 'std::vector<float>'):
        return t, 'series'
    if t in ('std::vector<int>', 'std::vector<std::size_t>', 'std::vector<unsigned>',
             'std::vector<std::uint64_t>'):
        return t, 'int_series'
    if t in ('std::vector<std::vector<double>>', 'Matrix'):
        return 'Matrix', 'matrix'
    if t == 'std::vector<std::complex<double>>':
        return 'ComplexField', 'complex_series'
    if t in ('std::string', 'std::vector<std::string>'):
        return t, 'text'
    if t.startswith('std::mt19937'):
        return t, 'seed'
    if t.startswith('std::function'):
        return t, 'callback'
    if t.startswith('std::vector<std::pair<'):
        return t, 'pair_series'
    if t in ('std::vector<std::vector<std::vector<double>>>',
             'std::vector<Matrix>'):
        return 'matrix_series', 'matrix_series'
    if t == 'std::vector<std::vector<std::complex<double>>>':
        return 'complex_matrix', 'complex_matrix'
    if t == 'std::vector<std::vector<std::pair<double, double>>>':
        return 'pair_matrix', 'pair_matrix'
    if t in enums:
        return t, 'enum'
    if t in structs:
        return t, 'struct'
    if t.startswith('std::vector<'):
        inner = t[t.index('<') + 1:t.rindex('>')]
        if inner in structs:
            return inner, 'struct_series'
    if t in ('json::Value', 'Value'):
        return t, 'text'
    raise Unsupported('%s (%s)' % (ptype, t))


def pointer_target(ptype):
    return '*' in ptype


def fmt_return(ret, expr, structs, enums, aliases):
    t = resolve(ret, aliases)
    if t in ('double', 'float', 'long double'):
        return 'fmt_number(%s)' % expr
    if t == 'bool':
        return '(%s ? std::string("true") : std::string("false"))' % expr
    if t in INT_TYPES:
        return 'fmt_integer(static_cast<double>(%s))' % expr
    if t == 'std::string':
        return 'fmt_text(%s)' % expr
    if t in ('std::vector<double>', 'std::vector<float>'):
        return 'fmt_series(%s)' % expr
    if t in ('std::vector<int>', 'std::vector<std::size_t>', 'std::vector<unsigned>',
             'std::vector<std::uint64_t>'):
        return 'fmt_int_series(%s)' % expr
    if t == 'Matrix':
        return 'fmt_matrix(%s)' % expr
    if t == 'ComplexField':
        return 'fmt_complex_series(%s)' % expr
    if t.startswith('std::pair<'):
        inner = t[t.index('<') + 1:t.rindex('>')]
        bits = [b.strip() for b in split_top(inner)]
        if len(bits) == 2 and all(b in ('double', 'float') for b in bits):
            return 'fmt_pair(%s)' % expr
        return 'fmt_label("%s")' % t
    if t == 'void':
        return None
    if t in enums:
        names = enums[t]
        return 'fmt_enum(static_cast<int>(%s), {%s})' % (
            expr, ', '.join('"%s"' % n for n in names))
    if t.startswith('std::vector<'):
        return 'fmt_label("%s")' % t
    members = structs.get(t, [])
    if not members:
        return 'fmt_label("%s")' % t
    lines = ['([&] {', '            std::string dump = "{";']
    for member in members:
        lines.append('            dump += "%s=";' % member['name'])
        lines.append('            dump += fmt_value(%s.%s);' % (expr, member['name']))
        lines.append('            dump += "; ";')
    lines.append('            dump += "}";')
    lines.append('            return dump;')
    lines.append('        }())')
    return '\n'.join(lines)


class Context:
    def __init__(self, structs, enums, aliases):
        self.structs = structs
        self.enums = enums
        self.aliases = aliases
        self.counter = 0
        self.prelude = []
        self.postlude = []
        self.series_slot = 0
        self.params_meta = []

    def next(self):
        value = self.counter
        self.counter += 1
        return value

    def number(self, ptype, name, default):
        slot = self.next()
        self.params_meta.append({'name': name, 'type': ptype, 'kind': 'number',
                                 'default': default})
        return slot

    def series_expr(self):
        slot = self.series_slot
        self.series_slot += 1
        if slot == 0:
            return 'in.time'
        if slot == 1:
            return 'in.series'
        return 'in.series'

    def build(self, param):
        ptype, kind = kind_of(param['type'], self.structs, self.enums, self.aliases)
        name = param['name']
        pointer = pointer_target(param['type'])

        if kind == 'number':
            slot = self.number(ptype, name, default_for(param))
            cast = SCALARS[ptype]
            expr = ('(in.numbers[%d] != 0.0)' % slot if cast == 'bool'
                    else 'static_cast<%s>(in.numbers[%d])' % (cast, slot))
            if pointer or ('&' in param['type'] and 'const' not in param['type']):
                self.prelude.append('        %s arg%d = %s;' % (ptype, slot, expr))
                self.postlude.append('        out += "  [out] %s = " + fmt_value(arg%d) + ";\\n";'
                                     % (name, slot))
                return ('&arg%d' % slot) if pointer else ('arg%d' % slot)
            return expr

        if kind == 'complex':
            slot = self.number(ptype, name, 1.0)
            self.params_meta[-1]['kind'] = 'complex'
            self.prelude.append(
                '        std::complex<double> arg%d(in.numbers[%d], in.numbers[%d]);'
                % (slot, slot, self.counter))
            self.counter += 1
            return 'arg%d' % slot

        if kind == 'series':
            expr = self.series_expr()
            if '&' in param['type'] and 'const' not in param['type']:
                slot = self.next()
                self.prelude.append('        std::vector<double> arg%d = %s;' % (slot, expr))
                self.postlude.append('        out += "  [out] %s = " + fmt_value(arg%d) + ";\\n";'
                                     % (name, slot))
                return 'arg%d' % slot
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'series', 'default': 0.0})
            return expr

        if kind == 'int_series':
            slot = self.next()
            element = 'int'
            if '<' in ptype:
                element = ptype[ptype.index('<') + 1:ptype.rindex('>')].strip()
            self.prelude.append('        std::vector<%s> arg%d = '
                                'to_index_series<%s>(in.series);' % (element, slot, element))
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'int_series', 'default': 0.0})
            return 'arg%d' % slot

        if kind == 'matrix':
            slot = self.next()
            name_tokens = set(param['name'].lower().split('_'))
            square = bool(name_tokens & {'weights', 'weight', 'covariance', 'correlation',
                                         'jacobian', 'transition', 'transitions', 'metric',
                                         'precision', 'kernel', 'graph', 'adjacency'})
            self.prelude.append('        Matrix arg%d = in.%s;'
                                % (slot, 'square' if square else 'matrix'))
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'matrix', 'default': 0.0})
            if '&' in param['type'] and 'const' not in param['type']:
                self.postlude.append('        out += "  [out] %s = " + fmt_value(arg%d) + ";\\n";'
                                     % (name, slot))
            return 'arg%d' % slot

        if kind == 'complex_series':
            slot = self.next()
            self.prelude.append('        std::vector<std::complex<double>> arg%d = '
                                'to_complex_series(in.series);' % slot)
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'complex_series', 'default': 0.0})
            return 'arg%d' % slot

        if kind == 'text':
            slot = self.next()
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'text', 'default': 0.0})
            if ptype in ('json::Value', 'Value'):
                self.prelude.append('        json::Value arg%d = json::Value::from_string(in.text);'
                                    % slot)
            elif ptype == 'std::vector<std::string>':
                self.prelude.append(
                    '        std::vector<std::string> arg%d = split_text(in.text);' % slot)
            else:
                self.prelude.append('        std::string arg%d = in.text;' % slot)
            return 'arg%d' % slot

        if kind == 'seed':
            slot = self.number(ptype, name, default_for(param))
            self.params_meta[-1]['kind'] = 'seed'
            self.prelude.append('        std::mt19937 arg%d(static_cast<unsigned>(std::fmax('
                                '1.0, in.numbers[%d])));' % (slot, slot))
            return 'arg%d' % slot

        if kind == 'enum':
            slot = self.number(ptype, name, default_for(param))
            self.params_meta[-1]['kind'] = 'enum'
            self.prelude.append('        auto arg%d = static_cast<climate::%s>(static_cast<int>('
                                'in.numbers[%d]));' % (slot, ptype, slot))
            return 'arg%d' % slot

        if kind == 'callback':
            slot = self.next()
            inner = ptype[ptype.index('<') + 1:ptype.rindex('>')]
            signature = re.match(r'^(.*?)\((.*)\)$', inner.strip())
            if signature:
                ret_type = signature.group(1).strip()
                arg_types = [t for t in split_top(signature.group(2)) if t != 'void']
            else:
                ret_type = inner.strip()
                arg_types = []
            lambda_params = ', '.join('%s a%d' % (t, i) for i, t in enumerate(arg_types))
            if ret_type == 'void':
                body = '{}'
            elif ret_type.startswith('std::pair'):
                body = '{ return std::make_pair(%s, %s); }' % (
                    'a0' if arg_types else '0.0', 'a1' if len(arg_types) > 1 else '0.0')
            elif ret_type == 'bool':
                body = ('{ return %s > 0.0; }' % (' + '.join('a%d' % i
                                                             for i in range(len(arg_types)))
                                                  or '0.0'))
            elif all(t in SCALARS for t in arg_types):
                body = ('{ return %s; }' % (' + '.join('a%d' % i
                                                        for i in range(len(arg_types)))
                                             or '0.0'))
            elif arg_types:
                body = '{ return static_cast<double>(a0.size()); }'
            else:
                body = '{ return 0.0; }'
            silenced = ''.join('(void)a%d; ' % i for i in range(len(arg_types)))
            body = body[:-1] + ' ' + silenced + '}' if body.startswith('{') else body
            self.prelude.append('        auto arg%d = [](%s) %s;' % (slot, lambda_params, body))
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'callback', 'default': 0.0})
            return 'arg%d' % slot

        if kind == 'pair_series':
            slot = self.next()
            element = 'double'
            inner = ptype[ptype.index('<') + 1:ptype.rindex('>')]
            parts = split_top(inner)
            if parts:
                bits = split_top(parts[0][parts[0].index('<') + 1:parts[0].rindex('>')])
                if len(bits) == 2:
                    element = bits[1]
            self.prelude.append('        std::vector<std::pair<%s, %s>> arg%d = '
                                'to_index_pairs<%s>(in.series);'
                                % (element, element, slot, element))
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'pair_series', 'default': 0.0})
            return 'arg%d' % slot

        if kind == 'complex_matrix':
            slot = self.next()
            self.prelude.append('        ComplexField arg%d;' % slot)
            self.prelude.append('        arg%d.push_back(to_complex_series(in.series));' % slot)
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'complex_matrix', 'default': 0.0})
            return 'arg%d' % slot

        if kind == 'pair_matrix':
            slot = self.next()
            self.prelude.append(
                '        std::vector<std::vector<std::pair<double, double>>> arg%d'
                '{to_pair_series(in.series)};' % slot)
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'pair_matrix', 'default': 0.0})
            return 'arg%d' % slot

        if kind == 'matrix_series':
            slot = self.next()
            self.prelude.append('        std::vector<Matrix> arg%d{in.matrix};' % slot)
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'matrix_series', 'default': 0.0})
            return 'arg%d' % slot

        if kind == 'struct':
            slot = self.next()
            self.prelude.append('        climate::%s arg%d{};' % (ptype, slot))
            self.fill_members('arg%d' % slot, ptype)
            if pointer:
                self.params_meta.append({'name': name, 'type': param['type'],
                                         'kind': 'struct', 'default': 0.0})
                return '&arg%d' % slot
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'struct', 'default': 0.0})
            if '&' in param['type'] and 'const' not in param['type']:
                self.postlude.append('        out += "  [out] %s = ";' % name)
                self.postlude.append('        out += fmt_value(arg%d);' % slot)
                self.postlude.append('        out += ";\\n";')
            return 'arg%d' % slot

        if kind == 'struct_series':
            slot = self.next()
            self.prelude.append('        std::vector<climate::%s> arg%d(3);' % (ptype, slot))
            self.prelude.append('        for (auto& element : arg%d) {' % slot)
            saved = self.prelude
            self.prelude = []
            self.fill_members('element', ptype, indent='            ')
            inner = self.prelude
            self.prelude = saved
            self.prelude.extend(inner)
            self.prelude.append('        }')
            self.params_meta.append({'name': name, 'type': param['type'],
                                     'kind': 'struct_series', 'default': 0.0})
            return 'arg%d' % slot

        raise Unsupported('%s -> %s' % (param['type'], kind))

    def fill_members(self, target, struct_name, indent='        '):
        members = self.structs.get(struct_name, [])
        used = 0
        for member in members:
            try:
                mtype, mkind = kind_of(member['type'], self.structs, self.enums, self.aliases)
            except Unsupported:
                continue
            if mkind == 'number':
                slot = self.number(mtype, member['name'], default_for({'name': member['name'],
                                                                       'default': None}))
                cast = SCALARS[mtype]
                value = ('(in.numbers[%d] != 0.0)' % slot if cast == 'bool'
                         else 'static_cast<%s>(in.numbers[%d])' % (cast, slot))
                self.prelude.append('%s%s.%s = %s;' % (indent, target, member['name'], value))
                used += 1
            elif mkind == 'series':
                self.prelude.append('%s%s.%s = in.series;' % (indent, target, member['name']))
                used += 1
            elif mkind == 'text':
                source = ('split_text(in.text)' if mtype == 'std::vector<std::string>'
                          else 'in.text')
                self.prelude.append('%s%s.%s = %s;' % (indent, target, member['name'], source))
                used += 1
        if used == 0 and indent == '        ':
            # Nothing scalar to fill: leave the default-constructed value.
            pass


def emit():
    functions, structs, enums, aliases, modules = scan()
    entries = []
    bodies = []
    used_ids = {}
    skipped = []
    for fn in functions:
        stage_id, stage_label = stage_of(fn['module'])
        formula_id = '%s::%s' % (fn['module'][:-4], fn['name'])
        if formula_id in used_ids:
            used_ids[formula_id] += 1
            formula_id = '%s#%d' % (formula_id, used_ids[formula_id])
        else:
            used_ids[formula_id] = 1
        ctx = Context(structs, enums, aliases)
        args = []
        try:
            for param in fn['params']:
                args.append(ctx.build(param))
        except Unsupported as error:
            skipped.append((formula_id, str(error), fn['ret']))
            continue

        qualifier = 'climate::'
        if fn['module'] == 'json.hpp':
            qualifier = 'climate::json::'
        call = '%s%s(%s)' % (qualifier, fn['name'], ', '.join(args))
        is_pointer = '*' in fn['ret']
        ret_fmt = fmt_return(fn['ret'], '(*result)' if is_pointer else 'result',
                             structs, enums, aliases)
        body = ['    case %d: {' % len(entries)]
        body.extend(ctx.prelude)
        if ret_fmt is None:
            body.append('        %s;' % call)
            body.append('        std::string out = "completed with no return value\\n";')
        elif is_pointer:
            body.append('        auto result = %s;' % call)
            body.append('        std::string out = result ? (%s) : std::string("null");'
                        % ret_fmt)
        else:
            body.append('        auto result = %s;' % call)
            body.append('        std::string out = %s;' % ret_fmt)
            body.append('        if (!out.empty() && out.back() != \'\\n\') out += "\\n";')
        body.extend(ctx.postlude)
        body.append('        return out;')
        body.append('    }')
        bodies.append(re.sub(r'in\.numbers\[(\d+)\]', r'num(in, \1)', '\n'.join(body)))

        entries.append({
            'id': formula_id,
            'name': fn['name'],
            'module': fn['module'],
            'stage': stage_id,
            'stage_label': stage_label,
            'signature': '%s %s(%s)' % (fn['ret'], fn['name'],
                                        ', '.join(p['raw'] for p in fn['params'])),
            'returns': fn['ret'],
            'params': ctx.params_meta,
        })

    if skipped:
        sys.stderr.write('skipped %d declarations:\n' % len(skipped))
        for item in skipped:
            sys.stderr.write('  %s :: %s -> %s\n' % item)

    out = []
    out.append('// Generated by tools/gen_formula_lab.py - do not edit by hand.')
    out.append('// Every public declaration in src/*.hpp has a call site below.')
    out.append('#include "formula_lab.hpp"')
    out.append('')
    for module in modules:
        out.append('#include "%s"' % module)
    out.append('')
    out.append('namespace climate {')
    out.append('namespace {')
    out.append('')
    out.append('double num(const LabInputs& in, std::size_t index) {')
    out.append('    return index < in.numbers.size() ? in.numbers[index] : 0.0;')
    out.append('}')
    out.append('')
    out.append('template <class I>')
    out.append('std::vector<I> to_index_series(const std::vector<double>& values) {')
    out.append('    std::vector<I> out;')
    out.append('    out.reserve(values.size());')
    out.append('    for (double value : values)')
    out.append('        out.push_back(static_cast<I>(std::llround(value)));')
    out.append('    return out;')
    out.append('}')
    out.append('')
    out.append('template <class I>')
    out.append('std::vector<std::pair<I, I>> to_index_pairs(const std::vector<double>& values) {')
    out.append('    std::vector<std::pair<I, I>> out;')
    out.append('    for (std::size_t i = 0; i + 1 < values.size(); ++i)')
    out.append('        out.emplace_back(static_cast<I>(std::llround(values[i])),')
    out.append('                     static_cast<I>(std::llround(values[i + 1])));')
    out.append('    return out;')
    out.append('}')
    out.append('')
    out.append('std::vector<std::pair<double, double>> to_pair_series('
               'const std::vector<double>& values) {')
    out.append('    std::vector<std::pair<double, double>> out;')
    out.append('    for (std::size_t i = 0; i + 1 < values.size(); ++i)')
    out.append('        out.emplace_back(values[i], values[i + 1]);')
    out.append('    return out;')
    out.append('}')
    out.append('')
    out.append('std::vector<std::complex<double>> to_complex_series('
               'const std::vector<double>& values) {')
    out.append('    std::vector<std::complex<double>> out;')
    out.append('    out.reserve(values.size());')
    out.append('    for (std::size_t i = 0; i < values.size(); ++i)')
    out.append('        out.emplace_back(values[i], 0.0);')
    out.append('    return out;')
    out.append('}')
    out.append('')
    out.append('std::vector<std::string> split_text(const std::string& text) {')
    out.append('    std::vector<std::string> out;')
    out.append('    std::string current;')
    out.append('    for (char c : text) {')
    out.append('        if (c == \',\') {')
    out.append('            out.push_back(current);')
    out.append('            current.clear();')
    out.append('        } else {')
    out.append('            current.push_back(c);')
    out.append('        }')
    out.append('    }')
    out.append('    if (!current.empty()) out.push_back(current);')
    out.append('    return out;')
    out.append('}')
    out.append('')
    out.append('}  // namespace')
    out.append('')
    out.append('const std::vector<LabFormula>& lab_formulas() {')
    out.append('    static const std::vector<LabFormula> table = [] {')
    out.append('        std::vector<LabFormula> list;')
    out.append('        list.reserve(%d);' % len(entries))
    out.append('        auto add = [&](std::size_t index, const char* id, const char* name,')
    out.append('                       const char* module, const char* stage,')
    out.append('                       const char* stage_label, const char* signature,')
    out.append('                       const char* returns) {')
    out.append('            LabFormula formula;')
    out.append('            formula.index = index;')
    out.append('            formula.id = id;')
    out.append('            formula.name = name;')
    out.append('            formula.module = module;')
    out.append('            formula.stage = stage;')
    out.append('            formula.stage_label = stage_label;')
    out.append('            formula.signature = signature;')
    out.append('            formula.returns = returns;')
    out.append('            list.push_back(formula);')
    out.append('        };')
    for position, entry in enumerate(entries):
        out.append('        add(%d, %s, %s, %s, %s, %s, %s, %s);' % (
            position, cstr(entry['id']), cstr(entry['name']), cstr(entry['module']),
            cstr(entry['stage']), cstr(entry['stage_label']),
            cstr(entry['signature']), cstr(entry['returns'])))
    out.append('        auto attach = [&](std::size_t index, std::size_t param_index, '
               'const char* name,')
    out.append('                           const char* type, const char* kind, double value) {')
    out.append('            LabParam param;')
    out.append('            param.name = name;')
    out.append('            param.type = type;')
    out.append('            param.kind = kind;')
    out.append('            param.default_value = value;')
    out.append('            list[index].params.resize(param_index + 1);')
    out.append('            list[index].params[param_index] = param;')
    out.append('        };')
    for position, entry in enumerate(entries):
        for index, param in enumerate(entry['params']):
            out.append('        attach(%d, %d, %s, %s, %s, %r);' % (
                position, index, cstr(param['name']), cstr(param['type']),
                cstr(param['kind']), param['default']))
    out.append('        return list;')
    out.append('    }();')
    out.append('    return table;')
    out.append('}')
    out.append('')
    out.append('std::string run_formula(std::size_t index, const LabInputs& in) {')
    out.append('    switch (index) {')
    out.extend(bodies)
    out.append('    default:')
    out.append('        return "no binding for this formula\\n";')
    out.append('    }')
    out.append('}')
    out.append('')
    out.append('}  // namespace climate')
    out.append('')
    open(os.path.join(ROOT, 'src', 'formula_lab_gen.cpp'), 'w', encoding='utf-8').write(
        '\n'.join(out))
    sys.stderr.write('generated %d runnable formulas (%d skipped)\n' % (
        len(entries), len(skipped)))
    return entries, skipped


def cstr(text):
    return '"%s"' % text.replace('\\', '\\\\').replace('"', '\\"')


if __name__ == '__main__':
    emit()
