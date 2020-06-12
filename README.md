# Adjust-Sigmoid-Parameters

Script realizado en C#, el cual recibe mediciones experimentales cargadas ya sea durante la operación o al finalizarla.
Posteriormente ajusta los valores experimentales a una curva sigmoidea mediante el método del gradiente decresciente, utilizando como función objetivo el mínimo error cuadrado y 
calcula la intersección entre las rectas correspondientes a:

- El punto de inflexión de la curva.
- El valor promedio de las 3 primeras mediciones.

Posteriormente a eso calcula el valor equivalente a Amilosa presente en almidón.

Este código fue utilizado en [este trabajo](https://revistas.unc.edu.ar/index.php/FCEFyN/article/view/16779/23397)

![ajuste](/imagenes/ajuste.png)
![convergencia](/imagenes/convergencia.png)
