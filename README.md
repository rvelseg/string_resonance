# string_resonance

Authors : R. Velasco-Segura and Pablo L. Rendón

Affiliations : Grupo de Acústica y Vibraciones, Centro de Ciencias
Aplicadas y Desarrollo Tecnológico (CCADET), Posgrado en Ciencias
Físicas (PCF), Universidad Nacional Autónoma de México (UNAM).

Source repository : https://github.com/rvelseg/string_resonance

## Description

A finite differences python script to simulate resonance in a string
with an obstacle. It forces the string at the left end with a sine,
which frequency is monotonically changing. This is, a linear or
exponential chirp.

Notice that the frequencies of resonance doesn't have wavelengths which
are perfect multiples of the double of the resonator length. Which
would be the case if the obstacle were not moving.

## Sample results

https://youtu.be/65gxCDMYnOc

## Acknowledgments

Thanks to Marcos Duarte for the peak detection script, which is
redistributed here under the terms of its MIT license,
https://github.com/demotu/BMC

This script started as a modified version of the code found in the
following forum http://stackoverflow.com/questions/10081942/ . Thanks
to anzupe, j08691, and kame.
