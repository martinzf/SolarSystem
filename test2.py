import datetime as dt
import numpy as np

today = dt.date.today()
old = dt.date(1800,3,21)

print((today - old).days)
print((old - today).days)
print(np.sign(1+2))