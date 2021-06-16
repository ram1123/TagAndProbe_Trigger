# -*- coding: utf-8 -*-
# @Author: Ram Krishna Sharma
# @Date:   2021-04-28
# @Last Modified by:   Ram Krishna Sharma
# @Last Modified time: 2021-05-23
import uproot

file = uproot.open('TnP_ntuple.root')

# print file["ntupler"]
# print file.keys()
# print file.classnames()
# print file["ntupler/EventTree"].keys()
# print file["ntupler"].classnames()

tree = file["ntupler/EventTree"]

print(tree)

print(tree.keys())

# print("|{0:10} | {1:10} | {2:10} | {3:40} | {4:15}".format("Event No.","Ele Pos.","Filter No.","Filter Name","Filter Decision"))
for i,name in enumerate(tree.arrays()['filterName32']):
    print("="*98)
    print("|{0:10} | {1:10} | {2:10} | {3:40} | {4:15}".format("Event No.","Ele Pos.","Filter No.","Filter Name","Filter Decision"))
    print("="*98)
    for j,Vecfilters in enumerate(name):
        for k,filters in enumerate(Vecfilters):
            print("|{0:10} | {1:10} | {2:10} | {3:40} | {4:15}".format(i, j, k, filters,tree.arrays()['filterDecision32'][i][j][k]))
