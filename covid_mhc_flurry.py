#!/usr/bin/env python3
"""
covid_mhc_flurry.py

runs mhc flurry on peptides

MIT License

Copyright (c) 2020 pdxgx

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

infile = "covid.pep"
peplist = [line.rstrip('\n') for line in open(infile)]
infile = mhcflurry_alleles.txt
allele_list = [line.rstrip('\n') for line in open(infile)]

with open("mhcflurry_input.csv", 'w', newline="") as myfile:
    combinations = [[y,x] for y in peplist for x in allele_list]
    wr = csv.writer(myfile)
    wr.writerows(combinations)