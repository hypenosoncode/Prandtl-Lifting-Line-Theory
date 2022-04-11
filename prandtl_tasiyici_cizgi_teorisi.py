########################################################################
"""
2021-2022 BAHAR DÖNEMİ UCK302 AERODİNAMİK-II ÖDEVİ 

Kodu Yazan: İsmail Selçuk ÇINAR

Tarih: 26/03/2022

12.GRUP ÜYELERİ:
190414151 - Gamze AYTEP
200414359 - İsmail Selçuk ÇINAR
200414361 - Arif Berk ARI
210414983 - Mehmet Selim AKPUNAR
"""
########################################################################
import math
import numpy as np
import matplotlib.pylab as plt

# KANAT DEĞERLERİ - WING VALUES
N = 11 # İstasyon Sayısı - Number of Stations
S = 20.67 # Kanat Alanı [m^2] - Wing Area
AR = 9.08 # Açıklık Oranı - Aspect Ratio
TR = 0.80 # Kök Oranı - Taper Ratio
alpha_angle_of_attack = 10.0 # Hücum Açısı [derece] - Twist Angle
alpha_infinity = 5.98 # Taşıma Eğrisi Eğimi [derece] - Lift Curve Slope
alpha_zero_lift = -1.250 # Sıfır Taşıma Açısı [derece] - Zero Lift Angle of Attack

# HESAPLAMALAR - CALCULATIONS
b = math.sqrt(AR * S) # Kanat Açıklığı [m] - Wing Span
MAC = S / b # Ortalama Veter Uzunluğu [m] - Mean Aerodynamic Chord
root_chord = (1.5 * (1 + TR) * MAC) / (1 + TR + TR ** 2)  # Kök Kanat Uzunluğu [m] - Root Chord
theta = np.linspace((math.pi / (2 * N)), (math.pi / 2), N) # π/2N den π/2 ye kadar olan bir dizi oluşturuldu. - A sequence from π/2N to π/2 is created.
k = (b / 2) * np.cos(theta)
c = root_chord * (1 - (1 - TR) * np.cos(theta))

# MÜ DEĞERİNİ HESAPLAMA ALANI - FIELD OF CALCULATION OF MU VALUE
mu = alpha_infinity * c / (4 * b)
mu_calculation_2 = mu * (np.array(alpha_angle_of_attack) - alpha_zero_lift)
mu_calculation_2_rad = np.deg2rad(mu_calculation_2)

# FORMÜL KISMI - FORMULA PART
theta_list = []
for i in range(1, 2 * N + 1, 2):
    theta_calculation = np.sin(i * theta) * (1 + (mu * i) / (np.sin(list(theta))))
    theta_list.append(theta_calculation)

theta_equation = np.asarray(theta_list)
theta_equation_transpoz = np.transpose(theta_equation)
theta_equation_full = np.linalg.inv(theta_equation_transpoz) # Formül dizi olacak şekilde güzel bir şekilde düzenlendi. - The formula is nicely arranged to be an array.

answer_list = np.matmul(theta_equation_full, mu_calculation_2_rad) # Hesaplamaya dahil edilebilmesi için yeni bir dizi oluşturuldu. - A new array has been created so that it can be included in the calculation.
my_area_calculation = np.divide((4 * b), c) # Alanın değişimine göre hesaplanan değerler. - Values calculated based on the change of the field.

general_iter = 0
for i in range(0, 10):
    iter = (np.sin((2*i+1) * theta)) * answer_list[i] * my_area_calculation
    general_iter = general_iter + iter

# Hesaplamaya rho dahil edilmiştir. - Rho is included in the calculation.
general_rho = 0 # Minimum sürükleme için rho sıfır kabul edilmelidir. - For minimum drag, rho should be assumed to be zero.
for i in range(1, 10):
    rho = (i+1 * (answer_list[i])**2 / (answer_list[0])**2)
    general_rho = general_rho + rho

CL_calculation = (math.pi * AR * answer_list[0]) # Gerçek hesaplanan CL değeri - Actual calculated CL value
CDi_induced_calculation = CL_calculation ** 2 / (math.pi * AR) * general_rho # Gerçek hesaplanan CL indüklenmiş değeri. - Actual calculated CL induced value.
CL_induced_graph = general_iter ** 2 / (math.pi * AR) * general_rho # Grafik 'e bastılacak değerler. - Values to print to chart

# HESAPLANAN DEĞERLERİ BASTIRMA ALANI - SUPPRESSING THE CALCULATED VALUES AREA
print("Tasima Katsayisi - Lift Coefficient (CL): ", CL_calculation)
print("Induklenmis Surukleme Katsayisi - Induced Drag Coefficient (CDi): ", CDi_induced_calculation)

# GRAFİK OLUŞTURMAK İÇİN DİZİ OLUŞTURMA ALANI - SERIES CREATION AREA TO CREATE A GRAPHIC    
CL_array = np.append(0, general_iter)
CDi_induced_array = np.append(0, CL_induced_graph)

# X-EKSENİ UZUNLUĞU - X-AXIS LENGTH
x_eksen = [b / 2]
for i in range(0, N):
    x_eksen.append(k[i])

# GRAFİK BASTIRMA ALANI - GRAPHIC CREATING AREA
fig1 = plt.figure("AERODYNAMICS HOMEWORK [12.GROUP] - CL")
plt.plot(x_eksen, CL_array, marker=".")
plt.title("NACA 65(2)-215 / CL")
plt.xlabel("Kanat Yarı-Açıklık Konumu - Wing Half-Opening Position (m)")
plt.ylabel("Taşıma Katsayısı - Lift Coefficient (CL)")
plt.grid()
plt.show()

fig2 = plt.figure("AERODYNAMICS HOMEWORK [12.GROUP] - CDi")
plt.plot(x_eksen, CDi_induced_array, marker=".")
plt.title("NACA 65(2)-215 / CDi")
plt.xlabel("Kanat Yarı-Açıklık Konumu - Wing Half-Opening Position (m)")
plt.ylabel("İndüklenmiş Sürükleme Katsayısı - Induced Drag Coefficient (CDi)")
plt.grid()
plt.show()
