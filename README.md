<header>
<h1>Solar System</h1>
<p>Python project for displaying the Solar System at a given date. </p>
</header> 

The project animates the orbits of the Solar System's planets from 
their current position to their position at a user specified date. 

The data was obtained from E.M. Standish (1992) 
https://www.rschr.de/PRPDF/aprx_pos_planets.pdf
and adapted in order to simplify the code. Degrees were converted to
radians. Certain quantites and their rates of change were combined into
others in order to only track the planets' Kepler elements (e, a, I, Ω, ω, M):
<li> Argument of perihelion = longitude of perihelion - longitude of 
ascending node. </li>
<li> Mean anomaly = mean longitude - longitude of perihelion + (correction
terms). </li>