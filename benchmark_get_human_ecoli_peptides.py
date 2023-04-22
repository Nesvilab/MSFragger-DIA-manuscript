#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.

import os

from pyopenms import ProteaseDigestion, AASequence

os.chdir(r"G:\Dropbox\papers_Fengchao\msfragger_dia\script\results\benchmark")

digestion = ProteaseDigestion()
digestion.setEnzyme("Trypsin/P")
digestion.setMissedCleavages(2)

human_peptides = set()
ecoli_peptides = set()

organism = 0  # 1 = human, 2 = ecoli
for line in open("2022-02-18-reviewed-UP000005640-UP000000625.fas").readlines():
    if line.startswith(">"):
        if "OS=Homo sapiens" in line:
            organism = 1
        elif "OS=Escherichia coli (strain K12)" in line:
            organism = 2
        else:
            print("There is a protein without organism information: " + line)
            exit(1)
        continue

    result = []
    digestion.digest(AASequence.fromString(line), result, 7, 50)
    if line.startswith("M"):
        result2 = []
        digestion.digest(AASequence.fromString(line[1:]), result2, 7, 50)
        result.extend(result2)

    if organism == 1:
        for peptide in result:
            human_peptides.add(peptide.toString())
    elif organism == 2:
        for peptide in result:
            ecoli_peptides.add(peptide.toString())

overlap = ecoli_peptides.intersection(human_peptides)

with open("human_ecoli_overlap_peptides.txt", "w") as f:
    for s in overlap:
        f.write(s + "\n")

print("Human peptides: " + str(len(human_peptides)))
print("E. coli peptides: " + str(len(ecoli_peptides)))
print("Overlapping peptides: " + str(len(overlap)))
print("Done!")

