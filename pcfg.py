import copy
from ptree import PTree, Node
from pprint import pprint
from more_itertools import chunked

def build_tree(pointers, initial_cell, words):
    ptr = pointers.get(initial_cell)
    if ptr is None:  
        return Node(initial_cell[2], words[initial_cell[0]:initial_cell[1]]) # Add word to node
    if len(ptr) == 1:
        children = build_tree(pointers, initial_cell[:2] + ptr, words)       # Recursively build a junction
        return Node(initial_cell[2], [children])
    if len(ptr) == 3:                                                        # etc.
        children = [build_tree(pointers, (initial_cell[0], ptr[0], ptr[1]), words),build_tree(pointers, (ptr[0], initial_cell[1], ptr[2]), words)]
        return Node(initial_cell[2], children)

class PRule(object):
    def __init__(self, variable, derivation, probability):
        self.variable = str(variable)
        self.derivation = tuple(derivation)
        self.probability = float(probability)

    def derivation_length(self):
        return len(self.derivation)

    def __repr__(self):
        compact_derivation = ' '.join(self.derivation)
        return self.variable + ' -> ' + compact_derivation + ' (' + str(self.probability) + ')'

    def __eq__(self, other):
        try:
            return self.variable == other.variable and self.derivation == other.derivation
        except Exception:
            return False


class PCFG(object):
    def __init__(self, start_variable='S', rules=None):
        if rules is None:
            self.rules = {}
        else:
            self.rules = copy.deepcopy(rules)  # A dictionary that maps an str object to a list of PRule objects
        self.start = start_variable  # Start symbol of the grammar

    def add_rule(self, rule):
        '''
        Adds a rule to dictionary of rules in grammar.
        '''
        if rule.variable not in self.rules:
            self.rules[rule.variable] = []
        self.rules[rule.variable].append(rule)

    def remove_rule(self, rule):
        '''
        Removes a rule from dictionary of rules in grammar.
        '''
        try:
            self.rules[rule.variable].remove(rule)
        except KeyError:
            pass
        except ValueError:
            pass

    @property
    def nonterminals(self):
        result = []

        for var, rules in self.rules.items():
            for rule in rules:
                if rule.derivation_length() == 2:
                    result.append(rule)

        return result

    @property
    def singletons(self):
        result = []

        for var, rules in self.rules.items():
            for rule in rules:
                if rule.derivation_length() == 1 and rule.derivation[0] in self.rules.keys():
                    result.append(rule)

        return result
    
    def all_that_derive(self, seq):
        result = []

        for var, rules in self.rules.items():
            for rule in rules:
                if isinstance(seq, str):
                    if rule.derivation == (seq,):
                        result.append(rule)

        return result 

    def to_near_cnf(self):
        '''
        Returns an equivalent near-CNF grammar.
        '''

        # Define a duplicate grammar, with an additional S0 start rule that derives S with probability 1        
        near_cnf = PCFG(start_variable='S0', rules={**self.rules, 'S0': [PRule('S0', ('S',), 1)]})

        #region Eliminate E rules
        # Epsilon rules are the rules whose derivation is a tuple consisting only of the empty string
        epsilon_rules = list(filter(lambda x: x.derivation == ('', ), [r for v in near_cnf.rules.values() for r in v]))
        eliminated = []

        # Until we have eliminated all the epsilon rules...
        while epsilon_rules != []:
            epsilon_rule = epsilon_rules.pop()  # Grab an epsilon rule
            eliminated.append(epsilon_rule)     # Add it to the list of eliminated rules
            near_cnf.remove_rule(epsilon_rule)  # Remove it from the result grammar
        
        # Now, we normalize rule probabilities...
            for rule in near_cnf.rules[epsilon_rule.variable]:                        # For all the rules that derive from this epsilon rule...
                rule.probability = rule.probability / (1 - epsilon_rule.probability)  # Reassign probability

            # RHS rules are the rules that can derive to the current epsilon rule
            rhs_rules = list(filter(lambda x: epsilon_rule.variable in x.derivation, [r for v in near_cnf.rules.values() for r in v]))

            for rhs_rule in rhs_rules:                          # For all these rules...
                derivation_as_list = list(rhs_rule.derivation)  # Convert the immutable tuple to a list
                
                alternatives = [[]]
                m = 0
                
                while derivation_as_list:                               # As long as we have not exhausted the entire derivation...
                    if epsilon_rule.variable not in derivation_as_list: # And given that the variable exists in the derivation...
                        break

                    index = derivation_as_list.index(epsilon_rule.variable) # Memorize the first index of the variable in the derivation
                    m += 1
                    new_alternatives = []

                    for alternative in alternatives: # For all possible derivations containing the epsilon rule...
                        new_alternatives.append(alternative + derivation_as_list[:index])     # One derivation has the epsilon rule variable
                        new_alternatives.append(alternative + derivation_as_list[:index + 1]) # The other does not
                    
                    derivation_as_list = derivation_as_list[index + 1:] # Progress to the next chunk of the derivation
                    alternatives = new_alternatives                     # Update list of alternatives with new ones 
                
                alternatives = {tuple(alternative + derivation_as_list) # Convert to set
                                    for alternative in alternatives 
                                        if tuple(alternative + derivation_as_list) != rhs_rule.derivation}

                for alternative in alternatives:  
                    k = len(rhs_rule.derivation) - len(alternative) # k = number of removed tokens
                    new_rule = PRule(rhs_rule.variable,             # Create updated rule with new probability and derivation
                                    alternative,        
                                    rhs_rule.probability * (epsilon_rule.probability ** k) * (1 - epsilon_rule.probability) ** (m - k) 
                                    ) 

                    if alternative not in map(lambda x: x.derivation, self.rules[rhs_rule.variable]): # Add novel rules
                        near_cnf.add_rule(new_rule)
                    else:
                        index = near_cnf.rules[rhs_rule.variable].index(new_rule)                     # Update existing rules
                        near_cnf.rules[rhs_rule.variable][index].probability += new_rule.probability
                    
                    # Create the new rule
                    near_cnf.add_rule(PRule(rhs_rule.variable, rhs_rule.derivation, rhs_rule.probability * (1 - epsilon_rule.probability) ** m))

                # Create a list of rules that have an epsilon derivation in their RHS
                epsilon_rules = list(filter(lambda x: x.derivation == ('', ) and x not in eliminated, [r for v in near_cnf.rules.values() for r in v]))
            #endregion

        #region Modify long rules
        rule_index = 1
        long_rules = list(filter(lambda x: x.derivation_length() > 2, [r for v in near_cnf.rules.values() for r in v])) # Filter rules with long derivation
        while long_rules: # As long as we have long rules...
            long_rule = long_rules.pop() # Pick one
            near_cnf.remove_rule(long_rule) # Delete it from the grammar
            pairs = list(chunked(long_rule.derivation, 2)) # Divide the rule into pairs
            pair = pairs[0]   # Pick first pair
            pairs = pairs[1:] # Update rest
            new_rule = PRule(long_rule.variable + str(rule_index), pair, 1) # Create new rule from pair
            rule_index += 1
            long_rule.derivation = tuple([new_rule.variable, *long_rule.derivation[2:]]) # Update rule
            if long_rule.derivation not in [r for v in near_cnf.rules.values() for r in v]: # Avoid dupes
                near_cnf.add_rule(long_rule)
            near_cnf.add_rule(new_rule)
            long_rules = list(filter(lambda x: x.derivation_length() > 2, [r for v in near_cnf.rules.values() for r in v])) # Regenerate list
        #endregion
        
        #region Fix rules containing terminals
        categories = set(map(lambda x: (x, ), near_cnf.rules.keys())) # Convert str keys into tuples
        terminals = set()

        rules_to_fix = []

        for var, rules in near_cnf.rules.items(): # Collect terminals
            for rule in rules:
                if rule.derivation not in categories:
                    if rule.derivation_length() == 1:
                        terminals.add(rule.derivation)
                    elif any(map(lambda x: (x, ) not in categories, rule.derivation)): # We need to fix any rule containing a terminal.
                        rules_to_fix.append(rule)
    
        while rules_to_fix: # As long as we have rules to fix...
            rule = rules_to_fix.pop() # Pick one
            near_cnf.remove_rule(rule) # Remove from grammar
            l, r = rule.derivation # Destructure the derivation

            new_rule_l, new_rule_r = PRule(l.capitalize(), tuple(), 1), PRule(r.capitalize(), tuple(), 1)

            if (l, ) in terminals:
                new_rule_l.derivation = (l,)
            
            if (r, ) in terminals:
                new_rule_l.derivation = (r,)
            new_rule = PRule(rule.variable, (new_rule_l.variable, new_rule_r.variable), rule.probability)
            
            if new_rule_r.derivation and new_rule_l.derivation:        # Add rules
                near_cnf.add_rule(new_rule_l)
                near_cnf.add_rule(new_rule_r)
                near_cnf.add_rule(new_rule)
            elif new_rule_r.derivation and not new_rule_l.derivation:
                near_cnf.add_rule(new_rule_r)
                near_cnf.add_rule(new_rule)
            elif new_rule_l.derivation and not new_rule_r.derivation:
                near_cnf.add_rule(new_rule_l)
                near_cnf.add_rule(new_rule)
            else:
                continue
            
        #endregion
        #region Cull duplicate rules in the grammar
        for var, group in near_cnf.rules.items():
            rules = []
            for rule in group:
                if rule not in rules:
                    rules.append(rule)
                elif rule.derivation_length() > 2:
                    pass
                    
                else:
                    continue
            near_cnf.rules[var] = rules
        #endregion

        return near_cnf

    def cky_parser(self, string):
        '''
        Parses the input string given the grammar, using the probabilistic CKY algorithm.
        If the string has been generated by the grammar - returns a most likely parse tree for the input string.
        Otherwise - returns None.
        The CFG is given in near-CNF.
        '''
        near = self.to_near_cnf()

        for rule in near.rules[self.start]:
            if rule.derivation == () and string == "":
                return PTree(Node(self.start, [string]), rule.probability)

        table = {}    # Represent table as an associative array
        pointers = {} # etc.

        words = string.split(" ")

        
        for j in range(1, len(words) + 1):
            for rule in self.all_that_derive(words[j-1]): # Add all keys that derive word to table
                table[(j-1, j, rule.variable)] = rule.probability
            for i in range(j-1, -1, -1): # Split string to different divisions
                for k in range(i+1, j):
                    for rule in near.nonterminals: 
                        # If a key is not in the table, return zero, as in the array impl
                        if table.get((i, j, rule.variable), 0) < rule.probability * table.get((i, k, rule.derivation[0]), 0) * table.get((k, j, rule.derivation[1]), 0):
                            table[(i, j, rule.variable)] =  rule.probability * table.get((i, k, rule.derivation[0])) * table.get((k, j, rule.derivation[1]))
                            # Add junction pointer
                            pointers[(i, j, rule.variable)] = (k, rule.derivation[0], rule.derivation[1])
                # Derive singleton rules
                relevant_singletons = [r for r in near.singletons if table.get((i, j, r.variable), 0) < r.probability * table.get((i, j, r.derivation[0]), 0)]
                while relevant_singletons:
                    rel_sin = relevant_singletons.pop()
                    table[(i, j, rel_sin.variable)] = rel_sin.probability * table[(i, j, rel_sin.derivation[0])]
                    pointers[(i, j, rel_sin.variable)] = rel_sin.derivation
                    relevant_singletons = [r for r in near.singletons if table.get((i, j, r.variable), 0) < r.probability * table.get((i, j, r.derivation[0]), 0)]

        if not table.get((0, len(words) , near.start)):
            print( "Not parsed by G.")
            return None
        
        if table.get((0, len(words) , near.start)) < 2**-50:
            print( "Not parsed by G.")
            return None

        return PTree(build_tree(pointers, (0, len(words), self.start), words), table[(0, len(words), self.start)])

    def is_valid_grammar(self):
        '''
        Validates that the grammar is legal (meaning - the probabilities of the rules for each variable sum to 1).
        '''
        return all(map(lambda x: 1 >= x >= 1-0.0001 , [sum(map(lambda x: x.probability, r)) for r in self.rules.values()]))
    
    def adjust_near_cnf_ptree(self, ptree, changes):
        '''
        THIS METHOD IS RELEVANT ONLY FOR THE BONUS QUSETION.
        Adjusts a PTree derived by a grammar converted to near-CNF, to the equivalent PTree of the original grammar.
        '''
        # TODO implement this method (for the bonus question)
        pass


class PCFGChange(object):
    NEW_START = 'new_start'
    EPSILON_RULE = 'epsilon_rule'
    AUXILIARY = 'auxiliary'

    def __init__(self, rule, change_type, info=None):
        '''
        THIS CLASS IS RELEVANT ONLY FOR THE BONUS QUSETION.
        Documents the specific change done on a PCFG.
        '''
        assert change_type in (PCFGChange.NEW_START, PCFGChange.EPSILON_RULE, PCFGChange.AUXILIARY)
        self.rule = rule
        self.change_type = change_type
        self.info = info
