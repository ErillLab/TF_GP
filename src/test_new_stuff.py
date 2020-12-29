#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 12:12:13 2020

@author: elia
"""


from objects.organism_factory import OrganismFactory
import json

def read_json_file(filename: str) -> dict:
    """Reads a JSON file and returns a dictionary with the content

    Args:
        filename: Name of the json file to read

    Returns:
        Dictionary with the json file info
    """

    with open(filename) as json_content:
        return json.load(json_content)


JSON_CONFIG_FILENAME = "config.json"
config = read_json_file(JSON_CONFIG_FILENAME)

configOrganism = config["organism"]
configOrganismFactory = config["organismFactory"]
configConnector = config["connector"]
configPssm = config["pssm"]



# ---------------------------------TEST----------------------------------------

organism_factory = OrganismFactory(
        configOrganism, configOrganismFactory, configConnector, configPssm
    )

#------------------------------------------------------------------------------
# Test get_organism
#------------------------------------------------------------------------------

parent_1 = organism_factory.get_organism()  # Parent 1

# Create another organism
parent_2 = organism_factory.get_organism()  # Parent 2


#------------------------------------------------------------------------------
# Test 'print' method
#------------------------------------------------------------------------------

print("\n")
parent_1.print()  # Print parent_1
print("\n")
parent_2.print()  # Print parent_2




#------------------------------------------------------------------------------
# Test 'count' methods
#------------------------------------------------------------------------------

# (parent_1)
# Number of recognizers
print('Number of recognizers:', parent_1.count_recognizers())
# Same thing as
print('Number of recognizers:', len(parent_1.recognizers))

# Number of connectors
print('Number of connectors:', parent_1.count_connectors())
# Same thing as
print('Number of connectors:', len(parent_1.connectors))


# (parent_2)
print('Number of recognizers:', parent_2.count_recognizers())
print('Number of connectors:', parent_2.count_connectors())


#------------------------------------------------------------------------------
# Test get_children
#------------------------------------------------------------------------------

# Recombination process
child_1, child_2 = organism_factory.get_children(parent_1, parent_2)

# Print children
print("\n")
child_1.print()  # Print child_1
print("\n")
child_2.print()  # Print child_2


#------------------------------------------------------------------------------
# Test 'export' method
#------------------------------------------------------------------------------

parent_1.export("test_export_parent1.txt")
child_1.export("test_export_child1.txt")

#------------------------------------------------------------------------------
# Test json import/export
#------------------------------------------------------------------------------

# Three random organisms
org_1 = organism_factory.get_organism()
org_2 = organism_factory.get_organism()
org_3 = organism_factory.get_organism()

print("\n")
org_1.print()  # Print org_1
print("\n")
org_2.print()  # Print org_2
print("\n")
org_3.print()  # Print org_3

# Export organisms
organism_factory.export_organisms([org_1, org_2, org_3], "three_random_org_json_test.json")
# Import organisms
org_list = organism_factory.import_organisms("three_random_org_json_test.json")

print("\n")
org_list[0].print()  # Print org_1
print("\n")
org_list[1].print()  # Print org_2
print("\n")
org_list[2].print()  # Print org_3


#------------------------------------------------------------------------------
# Test PSSM scores
#------------------------------------------------------------------------------

new_pssm = organism_factory.create_pssm(4)  # PSSM

new_pssm.print()

new_pssm.get_score("ccag")
new_pssm.get_score("ctgg")


#------------------------------------------------------------------------------
# Test mutate organisms
#------------------------------------------------------------------------------

'''
In the cinfig.json file: turn on one of the three types of mutation at a time.
For example, to test deletion:

"MUTATE_PROBABILITY_NODE_MUTATION": 0,
"MUTATE_PROBABILITY_DELETE_RECOGNIZER": 1,
"MUTATE_PROBABILITY_INSERT_RECOGNIZER": 0,
'''

org_1 = organism_factory.get_organism()

print("\nBefore mutation:")
org_1.print()
org_1.mutate(organism_factory)
print("\nAfter mutation:")
org_1.print()



# Example of output from intelligent deletion (with notes)
'''
Before mutation:
***** Organism 11 *****
tgag
 m: 16 s: 20
atAc
 m: 49 s: 14
gctg
 m: 28 s: 34
CcAt
 m: 15 s: 41   <--- (connectors to merge)
agAT   <-------------- DELETED RECOGNIZER
 m: 19 s: 21   <--- (connectors to merge)
aGtc

After mutation:
***** Organism 11 *****
tgag
 m: 16 s: 20
atAc
 m: 49 s: 14
gctg
 m: 28 s: 34
CcAt
 m: 34 s: 46.06517122512408   <------- ADJUSTED (merged) CONNECTOR *
aGtc



* mu = 34 = 15 + 19
  sigma = 46.065... = (41**2 + 21**2)**(1/2)
'''


# Example of output from intelligent insertion (with notes)
'''
Before mutation:
***** Organism 18 *****
agTa
 m: 27 s: 20   <------------------ NEW RECOGNIZER WILL BE INSERTED HERE
aaga
 m: 45 s: 41
aacC
 m: 43 s: 35
agat
 m: 50 s: 19
ttaG
 m: 46 s: 35
tata
 m: 28 s: 24
taCt

After mutation:
***** Organism 18 *****
agTa
 m: 14.294117647058824 s: 7.6923076923076925   <--- ADJUSTED CONNECTORS *
ggta  <------------------------------------------ INSERTED RECOGNIZER
 m: 12.705882352941178 s: 18.461538461538463   <--- ADJUSTED CONNECTORS *
aaga
 m: 45 s: 41
aacC
 m: 43 s: 35
agat
 m: 50 s: 19
ttaG
 m: 46 s: 35
tata
 m: 28 s: 24
taCt



* 14.294... + 12.705... = 27
  (7.692**2 + 18.461**2)**(1/2) = 20
  => the spacing between agTa and aaga is preserved
'''
















