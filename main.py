from ptree import PTree, Node
from pcfg import PCFG, PRule
from pprint import pprint
import regex as re

PARSELEX = re.compile('\s*"(?P<lex>[a-zA-Z].*)"\s*\[(?P<probability>[01].\d*)\]\s*')
PARSEDER = re.compile('\s*(?P<der>[a-zA-Z\-\s]+)\s*\[(?P<probability>[01].\d*)\]\s*')
def load_grammar(path_to_grammar):
    with open(path_to_grammar) as grammar_file: # LOAD GRAMMAR FROM FILE.
                                                # TODO: condense into strings? idk
        g = PCFG()
        for line in grammar_file:
            if line.startswith("#"):
                continue
            if line.startswith("start_variable"):
                g.start = line.split(" ")[-1]
            if "->" in line:
                label, rest = line.split(" -> ")
                opts = rest.split("|")
                # print(opts)
                lex = filter(lambda x: x is not None, map(lambda x: re.match(PARSELEX, x), opts))
                der = filter(lambda x: x is not None, map(lambda x: re.match(PARSEDER, x), opts))
                # print(*map(lambda x: x.groupdict(), lex))
                # print(*map(lambda x: x.groupdict(), der))
                for l in lex:
                    # print(l.groupdict())
                    g.add_rule(PRule(label, (l.group("lex").strip(), ), l.group("probability")))
                for d in der:
                    # print(d.groupdict())
                    g.add_rule(PRule(label, tuple(d.group("der").strip().split(" ")), d.group("probability")))
        return g.to_near_cnf()
    
if __name__ == "__main__":
    near_cnf = load_grammar("./grammar.txt")
    with open("./data.txt") as data:
        for line in data:
            print(near_cnf.cky_parser(line.strip()))