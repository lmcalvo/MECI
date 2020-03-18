larynx$GrupoEdad <- cut(larynx$age, breaks = c(40, 55, 75, 90))  #Creando nueva variable
#Estimación Kaplan Meier para la función de supervivencia sin tener en cuenta el sexo
mortalidad <- survfit(Surv(time, delta) ~ stage, data = larynx, type = "kaplan-meier")
summary(mortalidad)

ggsurvplot(fit = mortalidad, data = larynx, title = "Curva de Supervivencia", 
           xlab = "Tiempo", ylab = "Probabilidad de supervivencia")

ggsurvplot(mortalidad, fun = "cumhaz", xlab = "Tiempo", censor = T, 
           ylab = "Riesgo Acumulado", title = "Riesgo Acumulado")

print(mortalidad.km,print.rmean=T)
print(mortalidad.km_sexo,print.rmean=T)