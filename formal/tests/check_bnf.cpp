// Checks of morphe/formal/bnf.hpp and of Earley's diagnostics (earley.hpp).
//
// - The atoms: the terminals are the coarsest partition of the bytes refining the grammar's
//   classes and literal letters (counted on overlapping classes).
// - Named trees: an arithmetic grammar's trees evaluated give the value of 400 random
//   expressions, so the tree is the derivation; hidden rules are spliced; texts are the spans.
// - The regular fragment: 80 random grammars of literals, classes, groups, alternation, ?, * and +
//   accept exactly what the same expression accepts as an ECMAScript regular expression
//   (std::regex), on every word over {a, b, c} up to length 6.
// - Diagnostics, exactly: on 60 random finite languages enumerated completely, the furthest
//   viable prefix of every word up to length 5 is the longest prefix of some word of the language,
//   and the terminals expected there are exactly the letters that continue it.
// - Messages: where, after what, what is unexpected and what could follow, literals shown from
//   where the text left them; errors in a grammar's text located by line and column.
// - JSON's numbers (RFC 8259 §6) against their regular expression on 20 000 random strings; a
//   start other than the first rule; escapes; a list of 5000 items in a chart linear in its
//   length.

#include <morphe/check.hpp>
#include <morphe/formal/bnf.hpp>
#include <morphe/formal/earley.hpp>

#include <cstdint>
#include <functional>
#include <regex>
#include <set>
#include <string>
#include <vector>

namespace mf = morphe::formal;

namespace {

// A small deterministic generator (splitmix64).
struct Rng {
  std::uint64_t s;
  std::uint64_t next() {
    std::uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  std::size_t below(std::size_t n) { return static_cast<std::size_t>(next() % n); }
};

// A random expression over {a, b, c}, written both as BNF and as an ECMAScript regex; `finite`
// leaves out * and +, and `words` is then its language, enumerated.
struct Expr {
  std::string bnf, regex;
  std::set<std::string> words;
};

std::set<std::string> product(const std::set<std::string>& x, const std::set<std::string>& y) {
  std::set<std::string> out;
  for (const auto& u : x)
    for (const auto& v : y) out.insert(u + v);
  return out;
}

Expr random_expr(Rng& r, int depth, bool finite) {
  const std::size_t kinds = depth <= 0 ? 2 : finite ? 5 : 7;
  switch (r.below(kinds)) {
    case 0: {  // a literal of one or two letters
      std::string w;
      const std::size_t n = 1 + r.below(2);
      for (std::size_t i = 0; i < n; i++) w += static_cast<char>('a' + r.below(3));
      return {"\"" + w + "\"", w, {w}};
    }
    case 1: {  // a class
      // [^b] holds every byte but b, so it is left out of the finite languages, which are
      // enumerated over {a, b, c}.
      static const char* classes[] = {"[ab]", "[bc]", "[a-c]", "[^b]"};
      const std::size_t k = r.below(finite ? 3 : 4);
      std::set<std::string> ws;
      for (char c = 'a'; c <= 'c'; c++) {
        const bool in = k == 0 ? c != 'c' : k == 1 ? c != 'a' : k == 2 ? true : c != 'b';
        if (in) ws.insert(std::string(1, c));
      }
      return {classes[k], classes[k], ws};
    }
    case 2: {  // a sequence
      Expr x = random_expr(r, depth - 1, finite), y = random_expr(r, depth - 1, finite);
      return {"(" + x.bnf + " " + y.bnf + ")", "(?:" + x.regex + y.regex + ")", product(x.words, y.words)};
    }
    case 3: {  // an alternation
      Expr x = random_expr(r, depth - 1, finite), y = random_expr(r, depth - 1, finite);
      std::set<std::string> ws = x.words;
      ws.insert(y.words.begin(), y.words.end());
      return {"(" + x.bnf + " | " + y.bnf + ")", "(?:" + x.regex + "|" + y.regex + ")", ws};
    }
    case 4: {  // optional
      Expr x = random_expr(r, depth - 1, finite);
      std::set<std::string> ws = x.words;
      ws.insert("");
      return {"(" + x.bnf + ")?", "(?:" + x.regex + ")?", ws};
    }
    case 5: {
      Expr x = random_expr(r, depth - 1, finite);
      return {"(" + x.bnf + ")*", "(?:" + x.regex + ")*", {}};
    }
    default: {
      Expr x = random_expr(r, depth - 1, finite);
      return {"(" + x.bnf + ")+", "(?:" + x.regex + ")+", {}};
    }
  }
}

void all_words(std::size_t max_length, const std::function<void(const std::string&)>& f) {
  std::vector<std::string> level{""};
  f("");
  for (std::size_t n = 1; n <= max_length; n++) {
    std::vector<std::string> next;
    for (const auto& w : level)
      for (char c = 'a'; c <= 'c'; c++) next.push_back(w + c);
    for (const auto& w : next) f(w);
    level = std::move(next);
  }
}

// A random expression with its value.
std::pair<std::string, long long> random_arithmetic(Rng& r, int depth) {
  if (depth <= 0 || r.below(3) == 0) {
    const long long v = static_cast<long long>(r.below(10));
    return {std::to_string(v), v};
  }
  const auto [a, x] = random_arithmetic(r, depth - 1);
  const auto [b, y] = random_arithmetic(r, depth - 1);
  switch (r.below(3)) {
    case 0: return {"(" + a + "+" + b + ")", x + y};
    case 1: return {"(" + a + "-" + b + ")", x - y};
    default: return {"(" + a + ")*(" + b + ")", x * y};
  }
}

bool throws_at(const std::string& text, std::size_t line, std::size_t column) {
  try {
    mf::compile_bnf(text);
  } catch (const mf::BnfError& e) {
    return e.line == line && e.column == column;
  }
  return false;
}

}  // namespace

int main() {
  morphe::Suite suite;

  // ---- atoms --------------------------------------------------------------------------------------
  {
    const auto a = mf::compile_bnf("s ::= [0-9]+ (\"+\" [0-9]+)* ;");
    const auto b = mf::compile_bnf("s ::= [a-z]* [a-f0-9]? ;");
    suite.require("atoms: digits, '+' and the rest (3)", a.atoms() == 3 && a.atom('1') == a.atom('7') && a.atom('1') != a.atom('+'),
                  std::to_string(a.atoms()));
    suite.require("atoms of overlapping classes: a–f, g–z, 0–9 and the rest (4)",
                  b.atoms() == 4 && b.atom('a') == b.atom('f') && b.atom('f') != b.atom('g') && b.atom('0') == b.atom('9') && b.atom('0') != b.atom('a'),
                  std::to_string(b.atoms()));
    suite.require("a class that is one atom is one terminal", a.grammar().terminals == 3);
  }

  // ---- named trees --------------------------------------------------------------------------------
  {
    const auto g = mf::compile_bnf(R"G(
      # arithmetic, with the operators visible and the parentheses hidden
      expr   ::= term (op_add term)* ;
      op_add ::= op ;
      op     ::= "+" | "-" ;
      term   ::= factor (_times factor)* ;
      _times ::= mul ;
      mul    ::= "*" ;
      factor ::= number | "(" expr ")" ;
      number ::= [0-9]+ ;
    )G");
    const auto p = g.parse("12+3*(4-5)");
    bool shape = p.accepted && p.tree.name == "expr" && p.tree.children.size() == 3 && p.tree.children[0].text == "12" &&
                 p.tree.children[1].name == "op_add" && p.tree.children[2].name == "term" && p.tree.children[2].text == "3*(4-5)";
    suite.require("the named tree of 12+3*(4-5): expr over term, op_add, term", shape);
    // Evaluate with op_add's op and term's mul read as operators.
    const std::function<long long(const mf::ParseNode&)> eval = [&](const mf::ParseNode& n) -> long long {
      if (n.name == "number") return std::stoll(n.text);
      if (n.name == "factor") return n.child("number") ? eval(*n.child("number")) : eval(*n.child("expr"));
      long long v = 0;
      char op = '+';
      bool first = true;
      for (const auto& c : n.children) {
        if (c.name == "op_add") { op = c.text[0]; continue; }
        if (c.name == "mul") { op = '*'; continue; }
        const long long x = eval(c);
        if (first) { v = x; first = false; }
        else if (op == '+') v += x;
        else if (op == '-') v -= x;
        else v *= x;
      }
      return v;
    };
    suite.require("12+3*(4-5) evaluates to 9", p.accepted && eval(p.tree) == 9);
    const auto spliced = p.tree.find_all("mul");
    suite.require("the hidden rule _times is spliced: its mul is a child of term", spliced.size() == 1 && p.tree.children[2].child("mul") != nullptr);

    Rng r{7};
    bool all = true;
    std::size_t tested = 0;
    for (int i = 0; i < 400; i++) {
      const auto [text, value] = random_arithmetic(r, 4);
      const auto q = g.parse(text);
      all = all && q.accepted && eval(q.tree) == value && q.tree.text == text;
      tested++;
    }
    suite.require("400 random expressions: their trees evaluate to their values", all, std::to_string(tested));

    const auto t = g.parse("3*4", "term");
    suite.require("a start other than the first rule: 3*4 from term", t.accepted && t.tree.name == "term" && !g.parse("1+2", "term").accepted);
  }

  // ---- the regular fragment against std::regex ----------------------------------------------------
  {
    Rng r{2026};
    bool agree = true;
    std::size_t words = 0, accepted = 0;
    std::string first_disagreement;
    for (int k = 0; k < 80; k++) {
      const Expr e = random_expr(r, 3, false);
      const auto g = mf::compile_bnf("s ::= " + e.bnf + " ;");
      const std::regex re(e.regex, std::regex::ECMAScript);
      all_words(6, [&](const std::string& w) {
        const bool a = g.accepts(w), b = std::regex_match(w, re);
        if (a != b && first_disagreement.empty()) first_disagreement = e.bnf + " on '" + w + "'";
        agree = agree && a == b;
        words++;
        accepted += a;
      });
    }
    suite.require("80 random regular grammars accept what std::regex does, every word up to length 6", agree,
                  std::to_string(words) + " words, " + std::to_string(accepted) + " accepted" + (first_disagreement.empty() ? "" : "; " + first_disagreement));
  }

  // ---- diagnostics on finite languages, exactly ----------------------------------------------------
  {
    Rng r{11};
    bool furthest_ok = true, expected_ok = true;
    std::size_t tested = 0;
    for (int k = 0; k < 60; k++) {
      const Expr e = random_expr(r, 3, true);
      const auto g = mf::compile_bnf("s ::= " + e.bnf + " ;");
      const mf::Earley earley(g.grammar());
      std::set<std::string> prefixes;
      for (const auto& w : e.words)
        for (std::size_t i = 0; i <= w.size(); i++) prefixes.insert(w.substr(0, i));
      all_words(5, [&](const std::string& w) {
        std::size_t want = 0;
        for (std::size_t i = 0; i <= w.size(); i++)
          if (prefixes.count(w.substr(0, i))) want = i;
        std::vector<mf::Symbol> sym;
        for (char c : w) sym.push_back(g.atom(static_cast<unsigned char>(c)));
        const auto res = earley.parse(sym);
        furthest_ok = furthest_ok && res.furthest == want;
        std::set<mf::Symbol> next;
        const std::string head = w.substr(0, want);
        for (char c = 'a'; c <= 'c'; c++)
          if (prefixes.count(head + c)) next.insert(g.atom(static_cast<unsigned char>(c)));
        expected_ok = expected_ok && std::set<mf::Symbol>(res.expected.begin(), res.expected.end()) == next;
        tested++;
      });
    }
    suite.require("the furthest viable prefix is the longest prefix of a word of the language", furthest_ok, std::to_string(tested) + " words in 60 languages");
    suite.require("the terminals expected there are exactly the letters that continue it", expected_ok);
  }

  // ---- messages -----------------------------------------------------------------------------------
  {
    const auto g = mf::compile_bnf(R"(edit ::= "del" | "dup" | "ins" ;)");
    const auto x = g.parse("dux"), y = g.parse("d"), z = g.parse("");
    suite.require("'dux': at column 3, after \"du\", unexpected 'x', expected 'p' (of 'dup')",
                  !x.accepted && x.furthest == 2 && x.error == "at column 3, after \"du\": unexpected 'x'; expected 'p' (of 'dup')", x.error);
    suite.require("'d': ends early, expected 'el' (of 'del') or 'up' (of 'dup')",
                  !y.accepted && y.error == "the text ends after \"d\"; expected one of 'el' (of 'del') or 'up' (of 'dup')", y.error);
    suite.require("'': expected one of 'del', 'dup' or 'ins'", z.error == "the text is empty; expected one of 'del', 'dup' or 'ins'", z.error);
    const auto n = mf::compile_bnf("n ::= [0-9]+ (\".\" [0-9]+)? ;");
    const auto m = n.parse("12.x");
    suite.require("a class is shown as written: after \"12.\", expected [0-9]", m.error == "at column 4, after \"12.\": unexpected 'x'; expected [0-9]", m.error);
  }

  // ---- errors in a grammar's text ------------------------------------------------------------------
  {
    suite.require("a missing ';': the next rule's '::=' cannot continue the rule (line 2, column 3)", throws_at("a ::= \"x\"\nb ::= \"y\" ;", 2, 3));
    suite.require("an unknown name, at its place", throws_at("a ::= b ;", 1, 7));
    suite.require("a rule defined twice", throws_at("a ::= \"x\" ;\na ::= \"y\" ;", 2, 1));
    suite.require("an unterminated literal", throws_at("a ::= \"x ;", 1, 11));
    suite.require("an empty alternative", throws_at("a ::= \"x\" | ;", 1, 13));
    suite.require("a rule that derives no word", throws_at("a ::= \"x\" a ;", 1, 1));
    suite.require("an empty class", throws_at("a ::= [] ;", 1, 9));
  }

  // ---- JSON's numbers, escapes, long lists ----------------------------------------------------------
  {
    const auto g = mf::compile_bnf(R"(
      number   ::= "-"? int frac? exp? ;
      int      ::= "0" | [1-9] [0-9]* ;
      frac     ::= "." [0-9]+ ;
      exp      ::= [eE] [+\-]? [0-9]+ ;
    )");
    const std::regex re(R"(-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?)");
    Rng r{3};
    const std::string letters = "0123456789.-+eE";
    bool agree = true;
    std::size_t numbers = 0;
    for (int i = 0; i < 20000; i++) {
      std::string w;
      const std::size_t n = 1 + r.below(7);
      for (std::size_t k = 0; k < n; k++) w += letters[r.below(letters.size())];
      const bool a = g.accepts(w), b = std::regex_match(w, re);
      agree = agree && a == b;
      numbers += a;
    }
    suite.require("JSON's numbers: the grammar and the regular expression agree on 20 000 strings", agree, std::to_string(numbers) + " numbers");

    const auto e = mf::compile_bnf(R"(s ::= "\x41" [\]\-] "\\" "\"" ;)");
    suite.require("escapes: \\x41, \\] and \\- in a class, \\\\ and \\\"", e.accepts("A]\\\"") && e.accepts("A-\\\"") && !e.accepts("A^\\\""));

    const auto list = mf::compile_bnf(R"(list ::= item ("," item)* ; item ::= [a-z]+ ;)");
    std::string w = "a";
    for (int i = 1; i < 5000; i++) w += ",b";
    const mf::Earley earley(list.grammar());
    std::vector<mf::Symbol> sym;
    for (char c : w) sym.push_back(list.atom(static_cast<unsigned char>(c)));
    const auto res = earley.parse(sym);
    const auto p = list.parse(w);
    suite.require("a list of 5000 items: accepted, with 5000 item nodes", res.accepts && p.accepted && p.tree.all("item").size() == 5000);
    suite.require("its chart is linear in its length (fewer than 16 items a letter)", res.items < 16 * w.size(),
                  std::to_string(res.items) + " items for " + std::to_string(w.size()) + " letters");
  }

  return suite.done();
}
