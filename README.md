# eGFR operator

##### Description

`eGFR` operator calculates the glomerular filtration rate.

##### Usage

Input projection|.
---|---
`col` | 2nd factor age, 3rd factor gender
`y-axis, layer1`| is the value for serum creatine marker
`y-axis, layer2`| is the value for serum cystatin c marker


Output relations|.
---|---
`eGFR_0`| numeric, equation 0 calculation of the GFR estimate (e.g. per cell)
`eGFR_1`| numeric, equation 1 calculation of the GFR estimate (e.g. per cell)
`eGFR_2`| numeric, equation 2 calculation of the GFR estimate (e.g. per cell)
`eGFR_3`| numeric, equation 3 calculation of the GFR estimate (e.g. per cell)

##### Details

`gfr` operator estimates Glomerular filtration rate (GFR) is the best overall index of kidney function.
Four equations are used to calculate four eGFR values.

equation 0:
```r
    a <- ifelse (gender=="Female", -0.248, -0.207) 
    k <- ifelse (gender=="Female", 0.7, 0.9) 
    eGFR <-  135 * (min(crt/k, 1))^a * (max(crt/k, 1))^-0.601 * (min(cyt/0.8, 1))^-0.375 * (max(cyt/0.8, 1))^-0.711 * 0.995^age
    eGFR <- ifelse (gender=="Female",0.969 * eGFR, eGFR)
    eGFR <- ifelse (race=="Black",1.08 * eGFR, eGFR)
```

equation 1:
```r
    if (gender == "Female"){
           if   (crt <= 0.7) { eGFR <- 144 * (crt/ 0.7)^-0.329  * 0.995^age }
           if   (crt >  0.7) { eGFR <- 144 * (crt/ 0.7)^-1.209 * 0.995^age }
    }
    if (gender == "Male"){
           if   (crt <= 0.9) { eGFR <- 144 * (crt/ 0.9)^-0.411  * 0.995^age }
           if   (crt >  0.9) { eGFR <- 144 * (crt/ 0.9)^-1.209 * 0.995^age }
    }
    eGFR <- ifelse (race=="Black",1.159 * eGFR, eGFR)
```

equation 2:
```r
      if   (cyt <= 0.8) { eGFR <- 133 * (cyt/ 0.8)^-0.449 * 0.996^age }
      if   (cyt >  0.8) { eGFR <- 133 * (cyt/ 0.8)^-1.328 * 0.996^age }
      
    eGFR <- ifelse (gender=="Female",0.932 * eGFR, eGFR)
```

equation 3:
```r
    if (gender == "Female"){
      if   (crt <= 0.7 || cyt <= 0.8) { eGFR <- 130 * (crt/ 0.7)^-0.248 * (cyt/ 0.8)^-0.375 * 0.995^age }
      if   (crt <= 0.7 || cyt >  0.8) { eGFR <- 130 * (crt/ 0.7)^-0.248 * (cyt/ 0.8)^-0.711 * 0.995^age }
      if   (crt >  0.7 || cyt <= 0.8) { eGFR <- 130 * (crt/ 0.7)^-0.601 * (cyt/ 0.8)^-0.375 * 0.995^age }
      if   (crt >  0.7 || cyt >  0.8) { eGFR <- 130 * (crt/ 0.7)^-0.601 * (cyt/ 0.8)^-0.711 * 0.995^age }
    }
    if (gender == "Male"){
      if   (crt <= 0.9 || cyt <= 0.8) { eGFR <- 130 * (crt/ 0.7)^-0.207 * (cyt/ 0.8)^-0.375 * 0.995^age }
      if   (crt <= 0.9 || cyt >  0.8) { eGFR <- 130 * (crt/ 0.7)^-0.207 * (cyt/ 0.8)^-0.711 * 0.995^age }
      if   (crt >  0.9 || cyt <= 0.8) { eGFR <- 130 * (crt/ 0.7)^-0.601 * (cyt/ 0.8)^-0.375 * 0.995^age }
      if   (crt >  0.9 || cyt >  0.8) { eGFR <- 130 * (crt/ 0.7)^-0.601 * (cyt/ 0.8)^-0.711 * 0.995^age }
    }
      eGFR <- ifelse (race=="Black",1.159 * eGFR, eGFR)
```

##### References

ref for equation 0: [https://www.kidney.org/professionals/kdoqi/gfr_calculator]

