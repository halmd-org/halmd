Formulae and units
******************

Lennard-Jones
=============

Propagation laws
----------------

.. glossary::

 velocity-Verlet algorithm
  .. math::
    :nowrap:

    \begin{align*}
      \vec{r}_i(t+\Delta t) &= \vec{r}_i(t) + \vec{v}_i(t+\frac{\Delta t}{2})\Delta t \\
      \vec{v}_i(t+\frac{\Delta t}{2}) &= \vec{v}_i(t) +
      \frac{\vec{F}_i(t)}{m}\frac{\Delta t}{2} \\
      \vec{v}_i(t+\Delta t) &= \vec{v}_i(t+\frac{\Delta t}{2}) +
      \frac{\vec{F}_i(t+\Delta t)}{m}\frac{\Delta t}{2}
    \end{align*}


Interaction laws
----------------

.. glossary::

 Lennard-Jones potential
  .. math::
    U(r_{ij}) = 4\epsilon
    \Bigl(\bigl(\frac{\sigma}{r_{ij}}\bigr)^{12}
    - \bigl(\frac{\sigma}{r_{ij}}\bigr)^6\Bigr)

 shifted, truncated Lennard-Jones potential
  .. math::
    U_s(r_{ij}) =
    \begin{cases}
      4\epsilon\left(\left(\dfrac{\sigma}{r}\right)^{12} -
      \left(\dfrac{\sigma}{r}\right)^{6}\right)-U(r_c)\,, & r < r_c \,,\\
      0\,, & r > r_c \,.
    \end{cases}

 Lennard-Jones force
  .. math::
    \vec{F}(\vec{r}_{ij}) = \frac{48\epsilon}{\sigma}
    \Bigl(\bigl(\frac{\sigma}{r_{ij}}\bigr)^{13}
    - \frac{1}{2}\bigl(\frac{\sigma}{r_{ij}}\bigr)^7\Bigr) \hat{r}_{ij}


Base MD units
-------------

.. glossary::

 unit of length
  .. math:: [r] = \sigma

 unit of energy
  .. math:: [E] = \epsilon

 unit of time
  .. math:: [t] = \sqrt{\frac{m\sigma^2}{\epsilon}}

