import os

t = os.system("bigslice -i input --threshold 350 output350")
if t==1:
    t=os.system("bigslice -i input --threshold 450 output450")
if t==1:
    t=os.system("bigslice -i input --threshold 650 output650")
if t==1:
    t=os.system("bigslice -i input --threshold 750 output750")