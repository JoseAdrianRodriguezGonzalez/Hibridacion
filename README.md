# Hirbidacion
Proyecto de hibridación
>Este es el código del proyecto que gráfica la hibridación atómica.

>[!NOTE]
>Se encuentra en la versión 1.0.0, por lo que sigue en fase Beta

>[!WARNING]
>El código presenta errores en la parte de la gráfica. Se deben de realizar los siguiente cambios
>
> - En el cálculo de las coordendas cartesianas utilizar el radio atómica de Bohr, siendo de esta manera que:
> - 
>
>   $$x=\sin(\theta) *\cos(\phi) * ||r||$$
>
>   $$y=\sin(\theta) * \cos(\phi) * ||r||$$
>
>   $$z=\cos(\theta)*||r||$$
>> En el código se utilizó como magnitud el resultado de
>>
>> $$Y=|\theta,\phi|\^2$$
>>
>> Por lo que termina siendo un radio extremadamente grande.
>> - El siguiente error que existe es el uso de $\theta$ y de $\phi$, ya que solamente se está utilizando los valores pertenecientes de cada campo de forma general, sabiendo que:
>>
>>   $$\theta \epsilon [0,\pi] ^ \phi \epsilon [0,2*\pi]$$
>>   Por lo que cada calculo de coordenadas cartesianas en realidad, debería de tener el calculio para lo armónicos individuales por cada hibridación y no solo el rango del campo esférico
>>   
>>   ![Representación de las coordenadas esféricas](https://github.com/JoseAdrianRodriguezGonzalez/Hibridacion/blob/main/assets/The-point-A-on-spherical-plane_Q320.jpg)
>>   
>> - El tercer error es de que en la realización de cada armónico esférico, no se considera que es una magintud, es decir
>>   
>>  $$Y|\theta,\phi|\^2$$
>>    Por lo que se debe corregir la falta del cuadrado

>[!TIP]
>Al terminar de corregir estos errores se puede realizar una interfaz que sea amigable con los usuarios
