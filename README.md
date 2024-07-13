# Hibridación
El proyecto está seccionado en diferentes etapas para gráficar la hibridación
 * [La parte radial](https://github.com/JoseAdrianRodriguezGonzalez/Hibridacion/tree/2bf4264a83fbd573a454ea5dacbccdfb4fb00c84/radial)
 * [Los armónicos esféricos reales](https://github.com/JoseAdrianRodriguezGonzalez/Hibridacion/tree/2bf4264a83fbd573a454ea5dacbccdfb4fb00c84/Real)
 * [Los armónicos esféricos complejos](https://github.com/JoseAdrianRodriguezGonzalez/Hibridacion/tree/2bf4264a83fbd573a454ea5dacbccdfb4fb00c84/spherical_complex)
 * [Funciones de probabilidad](https://github.com/JoseAdrianRodriguezGonzalez/Hibridacion/tree/2bf4264a83fbd573a454ea5dacbccdfb4fb00c84/wavefunction_Probability)
   - [2D](https://github.com/JoseAdrianRodriguezGonzalez/Hibridacion/tree/2bf4264a83fbd573a454ea5dacbccdfb4fb00c84/wavefunction_Probability/2D)
   - [3D](https://github.com/JoseAdrianRodriguezGonzalez/Hibridacion/tree/2bf4264a83fbd573a454ea5dacbccdfb4fb00c84/wavefunction_Probability/3D)
 * [Hibridiación.](https://github.com/JoseAdrianRodriguezGonzalez/Hibridacion/tree/2bf4264a83fbd573a454ea5dacbccdfb4fb00c84/hibridacion)

Cada tema contiene una pequeñas explicación de la matemática empleada
# Instalación
Para poder instalar el proyecto necesita contar con la versión más reciente de jupyter, o por lo menos la versión 6 de jupyter, en dado caso puede visualizarlo en Google Colab, sin embargo, tomese en cuenta el *warning* de abajo.

>[!NOTE]
>En está sección principal, contendrá una breve explicación del proyecto a su vez que contendrá imágenes

>[!WARNING]
>La sección de hibridación, no se puede visualizar por completo en Colab, por lo que se recomienda que se pueda visualizar en algún software que permita visualizar archivos de Jupyter

# Colaboradores

Para hacer posible este proyecto fue realizado por:

* [Isabel](https://github.com/Isa-Rocha31415)
* [Adrián](https://github.com/JoseAdrianRodriguezGonzalez)

 ## Ecuación de Erwin Schrödinger
 Las gráficas provistas utilizan la solución de la ecuación de Schrödinger para el átomo del hidrógeno.

 $$H\psi=i\hbar\frac{\partial \psi}{\partial t}$$

>>   ![Armónicos esféricos](https://github.com/JoseAdrianRodriguezGonzalez/Hibridacion/blob/main/assets/A4.png)
  
## Solución de la ecuación de Erwin Schrödinger
Siendo:
$$i\hbar \frac{\partial \Psi}{\partial t} = -\frac{\hbar^2}{2m}\nabla^2 \Psi+V(x,y,z)\Psi$$
Podemos observar que la función $\Psi(x,y,z,t)\rightarrow \psi(x,y,z)T(t)$
$\therefore$
$$i\hbar \frac{\partial \psi T}{\partial t} = -T\frac{\hbar^2}{2m}\nabla^2 \psi+V(x,y,z)\psi T$$