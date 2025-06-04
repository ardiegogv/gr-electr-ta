# Caracterización de un Grupo Electrógeno

**Caracterización energética y ambiental de un grupo electrógeno de 6kVA**

Este repositorio contiene el código desarrollado para el análisis experimental de un grupo electrógeno en el contexto del laboratorio de Termodinámica Aplicada. El script realiza el procesamiento de datos, cálculos termodinámicos y estimaciones de impacto económico y ambiental bajo distintas condiciones de operación.

---

## Funcionalidades principales

- Cálculo del flujo másico de aire usando la norma de ISO 5167-2 para placa orificio.
- Determinación del flujo másico de diésel y potencia eléctrica generada.
- Evaluación del exceso de aire y composición química de gases de escape.
- Estimación de la temperatura adiabática de combustión utilizando CoolProp.
- Balance energético mediante:
  - Energías asociadas a reactantes y productos.
  - Diferencias de entalpía con propiedades termodinámicas.
- Comparación mensual entre el uso del grupo electrógeno y el Sistema Eléctrico Nacional (SEN):
  - Emisiones de CO₂.
  - Consumo y costos mensuales.

---

## Archivos

- `gr_electrógeno.py`: Script principal con todos los cálculos.

---

## Requisitos

- Librerías necesarias:
  - `numpy`
  - `CoolProp`

### Instalación rápida

```bash
pip install numpy CoolProp
```

---

## Datos de entrada por ensayo

Cada ensayo experimental se define como un diccionario con las siguientes variables:

- Volumen de diésel consumido (mL)  
- Tiempo de combustión (s)  
- Número de impulsos eléctricos generados  
- Tiempo de medición de impulsos (s)  
- Altura diferencial manométrica (cm H₂O)  
- Temperatura de gases de escape (°C)

---

## Resultados calculados por ensayo

- Flujo másico de aire y diésel.  
- Relación aire-combustible (RAC) y exceso de aire.  
- Composición estequiométrica y real de productos.  
- Composición seca de gases (como las mide un analizador).  
- Temperatura adiabática de combustión.  
- Potencia al eje, potencia eléctrica y rendimiento global.  
- Comparación entre métodos energéticos.

---

## Comparación mensual de operación

Se evalúa una demanda mensual de 3000 kWh para calcular:

- Litros de diésel requeridos  
- Emisiones mensuales de CO₂  
- Costos del combustible  
- Comparación con el SEN (emisiones y costos)

---

## Autor

**Diego A. García**  
Ingeniería Civil Aeroespacial – Universidad de Concepción  
[ar.diego.gv@gmail.com](mailto:ar.diego.gv@gmail.com)