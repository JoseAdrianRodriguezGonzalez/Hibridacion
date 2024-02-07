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
>> $$Y=|\theta,\phi|\^2
>>
>> Por lo que termina siendo un radio extremadamente grande.
>> - El siguiente error que existe es el uso de $\theta$ y de $\phi$
