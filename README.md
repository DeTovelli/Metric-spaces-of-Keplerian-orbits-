# Metric-spaces-of-Keplerian-orbits-


The repository contains Pascal-language software tools for calculating metric distances in Keplerian orbital spaces. The algorithms are based on theoretical methods of celestial mechanics and the work of K. V. Kholshevnikov.

The program compares the orbital elements of two celestial bodies at identical epochs and calculates five fundamental metrics of orbital proximity. | Программа сопоставляет орбитальные элементы двух небесных тел на одинаковые эпохи и рассчитывает пять фундаментальных метрик близости орбит.


## 🚀 Key Features | Ключевые возможности

* **Time Synchronization:** Automatically searches and compares orbital data on the same Julian Date (`JD`). | **Синхронизация по времени:** Автоматический поиск и сопоставление орбитальных данных на одинаковые Юлианские даты (`JD`).
* **Kholshevnikov metrics:** Calculation of distances in factor spaces of orbits (\(\rho_1\), \(\rho_2\), \(\rho_5\)) | **Метрики Холшевникова:** Расчет расстояний в фактор-пространствах орбит (\(\rho_1\), \(\rho_2\), \(\rho_5\)).
* **D-criteria:** Calculation of classical astronomical orbital similarity criteria (Southworth-Hawkins \(D_{SH}\) and Drummond \(D_D\) | **D-критерии:** Вычисление классических астрономических критериев сходства орбит (Саутворта-Хокинса \(D_{SH}\) и Драммонда \(D_D\)).

### Algorithm for calculating distances between orbits | Алгоритм расчета расстояний между орбитами

The program compares the orbital elements of two bodies at the same epochs (`JD1 == JD2`) and calculates the following metrics: | Программа сопоставляет орбитальные элементы двух тел на одинаковые эпохи (`JD1 == JD2`) и рассчитывает следующие метрики:

1. **Mutual inclination of orbits | Взаимное наклонение орбит ($I$):**
   
   $$\cos I = \cos i_1 \cos i_2 + \sin i_1 \sin i_2 \cos(\Omega_1 - \Omega_2)$$

3. **Orbital metrics of Kholshevnikov | Орбитальные метрики Холшевникова ($\rho_1, \rho_2, \rho_5$):**
   
   $$\rho_2 = \sqrt{(1+e_1^2)p_1 + (1+e_2^2)p_2 - 2\sqrt{p_1 p_2}(\cos I + e_1 e_2 \cos P)}$$

5. **THE SOUTHWORTH & HAWKINS D-CRITERION | D-критерий Саутворта-Хокинса ($D_{SH}$):**
   $$D_{SH} = \sqrt{(q_1 - q_2)^2 + (e_1 - e_2)^2 + 4\sin^2\frac{I}{2} + (e_1 + e_2)^2 \sin^2\Pi}$$

6. **D-критерий Драммонда ($D_D$):**
   $$D_D = \sqrt{\left(\frac{e_1 - e_2}{e_1 + e_2}\right)^2 + \left(\frac{q_1 - q_2}{q_1 + q_2}\right)^2 + \left(\frac{I}{\pi}\right)^2 + \left(\frac{e_1 + e_2}{2}\right)^2 \left(\frac{\theta}{\pi}\right)^2}$$

   ## 📊 Структура данных

### Входные файлы (`infile1` / `infile2`)
Текстовые файлы должны содержать строки со следующими параметрами, разделенными пробелами:
`JD` `t` `a` `e` `i` `om` `w` `M` `q`

Где:
* **JD** — Юлианская дата (используется для синхронизации потоков)
* **t** — Время / эпоха
* **a** — Большая полуось
* **e** — Эксцентриситет
* **i** — Наклонение орбит (градусы)
* **om** (\(\Omega\)) — Долгота восходящего узла (градусы)
* **w** (\(\omega\)) — Аргумент перицентра (градусы)
* **M** — Средняя аномалия (градусы)
* **q** — Перигелийное расстояние

### Выходной файл (`outfile`)
Результаты записываются в структурированном виде:
```text
[JD1] [t1] [ro1] [ro2] [ro5] [D_SH] [D_D]
```

---

## 🛠 Особенности реализации кода

В алгоритме применены важные оптимизации для повышения точности астрономических расчетов:
* **Безопасное сравнение дат:** Поиск совпадений `JD1` и `JD2` реализован через дельту (`Abs(JD1 - JD2) < 1e-6`), что предотвращает ошибки округления вещественных чисел (`Real`/`Double`).
* **Оптимизация прохода файлов:** При нахождении совпадения внутренний цикл немедленно прерывается (`Break`), переводя указатель на обработку следующей эпохи.
* **Защита от математических исключений (Clamping):** Для экстремальных орбит (с высоким эксцентриситетом $e \to 1$) тригонометрические функции защищены от ошибок округления чисел с плавающей точкой (`NaN`/`Domain Error`). Аргументы функций `ArcCos` и `ArcSin` принудительно ограничиваются диапазоном `[-1, 1]`, а также добавлена проверка деления на ноль при вычислении секанса взаимного наклонения.

## 💻 Требования и запуск

Для компиляции кода требуется любой современный компилятор Pascal:
* **Free Pascal Compiler (FPC)** 3.0+
* **Delphi** или **Lazarus IDE**

```bash
# Пример компиляции через FPC
fpc main.pas
./main
```


## 📄 Лицензия
Проект распространяется под лицензией MIT. Подробнее см. файл `LICENSE`
