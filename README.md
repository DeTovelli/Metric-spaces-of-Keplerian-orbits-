# Metric-spaces-of-Keplerian-orbits-
### Algorithm for calculating distances between orbits | Алгоритм расчета расстояний между орбитами

The program compares the orbital elements of two bodies at the same epochs (`JD1 == JD2`) and calculates the following metrics: | Программа сопоставляет орбитальные элементы двух тел на одинаковые эпохи (`JD1 == JD2`) и рассчитывает следующие метрики:

1. **Mutual inclination of orbits | Взаимное наклонение орбит ($I$):**
   
   $$\cos I = \cos i_1 \cos i_2 + \sin i_1 \sin i_2 \cos(\Omega_1 - \Omega_2)$$

3. **Orbital metrics of Kholshevnikov | Орбитальные метрики Холшевникова ($\rho_1, \rho_2, \rho_5$):**
   
   $$\rho_2 = \sqrt{(1+e_1^2)p_1 + (1+e_2^2)p_2 - 2\sqrt{p_1 p_2}(\cos I + e_1 e_2 \cos P)}$$

5. **THE SOUTHWORTH & HAWKINS D-CRITERION | D-критерий Саутворта-Хокинса ($D_{SH}$):**
   $$D_{SH} = \sqrt{(q_1 - q_2)^2 + (e_1 - e_2)^2 + 4\sin^2\frac{I}{2} + (e_1 + e_2)^2 \sin^2\Pi}$$

6. **D-критерий Драммонда ($D_D$):**
   $$D_D = \sqrt{\left(\frac{e_1 - e_2}{e_1 + e_2}\right)^2 + \left(\frac{q_1 - q_2}{q_1 + q_2}\right)^2 + \left(\frac{I}{\pi}\right)^2 + \left(\frac{e_1 + e_2}{2}\right)^2 \left(\frac{\theta}{\pi}\right)^2}$$
