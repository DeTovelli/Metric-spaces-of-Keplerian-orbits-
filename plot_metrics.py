pythonimport pandas as pd
import matplotlib.pyplot as plt
import os

# Проверяем наличие файла с результатами
if not os.path.exists('file.TXT'):
    print("Ошибка: Файл file.TXT не найден! Сначала запустите расчет на Pascal.")
    exit()

# Загружаем данные (разделитель — пробелы, пропускаем пустые строки)
data = pd.read_csv('file.TXT', sep=r'\s+', index_col=False)

# Настройка стиля графиков
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# 1. График метрик Холшевникова
ax1.plot(data['t'], data['ro1'], label=r'$\rho_1$', color='#1f77b4', alpha=0.8)
ax1.plot(data['t'], data['ro2'], label=r'$\rho_2$', color='#ff7f0e', alpha=0.8)
ax1.plot(data['t'], data['ro5'], label=r'$\rho_5$', color='#2ca02c', alpha=0.8)
ax1.set_title('Эволюция метрик Холшевникова во времени', fontsize=12, fontweight='bold')
ax1.set_ylabel('Расстояние между орбитами', fontsize=10)
ax1.legend(loc='upper right')

# 2. График D-критериев (Саутворта-Хокинса и Драммонда)
ax2.plot(data['t'], data['D_SH'], label='$D_{SH}$ (Southworth-Hawkins)', color='#d62728', linewidth=2)
ax2.plot(data['t'], data['D_D'], label='$D_D$ (Drummond)', color='#9467bd', linewidth=2)
ax2.set_title('Сравнение D-критериев сходства орбит', fontsize=12, fontweight='bold')
ax2.set_xlabel('Время / Эпоха (t)', fontsize=10)
ax2.set_ylabel('Значение критерия', fontsize=10)
ax2.legend(loc='upper right')

# Оптимизация полей и сохранение
plt.tight_layout()
output_image = 'orbit_metrics_plot.png'
plt.savefig(output_image, dpi=300)
print(f"Успех! График успешно сохранен как {output_image}")