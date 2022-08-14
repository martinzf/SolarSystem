# Solar System

### About
Python project for displaying the Solar System.
It animates the Solar System's planets from their current position to 
one at a user specified date between 3000 BC and 3000 AD.

### How to run
1. Either fork or download the project and open the folder in the 
command line.
2. Install dependencies using the `pip install -r requirements.txt` command.
3. Run the script with `python main.py`, or `python3 main.py` if both Python 2
and Python 3 are installed on your computer.
4. You will be prompted to input a date in ISO format, and whether it
is AD or BC. Do so and wait while the programme runs &mdash; it could take a few 
seconds.
5. A GIF showing the Solar System evolving in time from today to your chosen date
will be generated in the project directory with the name `solar_system.gif`
and opened.

### Preview
![alt text](preview.gif)

### References
The data was obtained from 
[E.M. Standish (1992)](https://www.rschr.de/PRPDF/aprx_pos_planets.pdf)
and adapted in order to simplify the code. Degrees were converted to
radians. Certain quantities and their rates of change were combined into
others in order to only track the planets' classical Keplerian orbital 
elements (e, a, ω, i, Ω, M):
- Argument of perihelion = longitude of perihelion - longitude of 
ascending node 
- Mean anomaly = mean longitude - longitude of perihelion + (correction
terms) 