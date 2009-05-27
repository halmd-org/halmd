Formulae and units
******************

Lennard-Jones
=============

Interaction laws
----------------

.. glossary::

 Lennard-Jones potential
  .. math::
    U(r_{ij}) = 4\epsilon
    \Bigl(\bigl(\frac{\sigma}{r_{ij}}\bigr)^{12}
    - \bigl(\frac{\sigma}{r_{ij}}\bigr)^6\Bigr)

 Shifted, truncated Lennard-Jones potential
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
    - \frac{1}{2}\bigl(\frac{\sigma}{r_{ij}}\bigr)^8\Bigr) \hat{r}_{ij}


Base MD units
-------------

.. glossary::

 length
  .. math:: [r] = \sigma

 energy
  .. math:: [E] = \epsilon

 time
  .. math:: [t] = \sqrt{\frac{m\sigma^2}{\epsilon}}

