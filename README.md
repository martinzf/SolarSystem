# Solar System

### About
Python 3.11 project with the purpose of practising animation.
It animates the Solar System's planets from their current position to 
one at a user specified date between 3000 BC and 3000 AD.

### How to run
1. Clone the project and open the folder from the 
command line.
1. (Optional) Create a virtual environment to install dependencies with the `python -m venv venv` command, followed by `venv/Scripts/activate`. (Replace `python` with `python3` in all commands if both Python 2 and Python 3 are installed on your computer).
1. Install dependencies using the `pip install -r requirements.txt` command.
1. Run the programme with `python main.py`.
1. You will be prompted to input a date in ISO format, and whether it
is AD or BC. Do so and wait while the programme runs &mdash; it could take a few 
seconds.
1. A GIF showing the Solar System evolving in time from today to your chosen date
will be generated in the project directory with the name `solar_system.gif`
and opened.

### Preview
![alt text](preview.gif)

### Reference
The astronomical data was obtained from 
[NASA's JPL](https://ssd.jpl.nasa.gov/planets/approx_pos.html), courtesy of E.M. Standish (1992). 