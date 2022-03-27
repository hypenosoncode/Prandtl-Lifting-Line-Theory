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

# DEĞERLER
N = 11 # İstasyon Sayısı
S = 20.67 # Kanat Alanı [m^2] - Wing Area
AR = 9.08 # Açıklık Oranı - Aspect Ratio
TR = 0.80 # Taper Ratio
alpha_hucum = 10.0 # Hücum Açısı [derece] - Twist Angle
alpha_sonsuz = 5.98 # Taşıma Eğrisi Eğimi [derece] - Lift Curve Slope
alpha_zero_lift = -1.250 # Sıfır Taşıma Açısı [derece] - Zero Lift Angle of Attack

# HESAPLAMALAR
b = math.sqrt(AR * S) # Kanat Açıklığı (m)
MAC = S / b # Ortalama Veter Uzunluğu (m)
kok_chord = (1.5 * (1 + TR) * MAC) / (1 + TR + TR ** 2)  # Kök kanat[m] - root chord

theta = np.linspace((math.pi / (2 * N)), (math.pi / 2), N) # Dizi oluşturuldu.

k = (b / 2) * np.cos(theta)
c = kok_chord * (1 - (1 - TR) * np.cos(theta))

mu = alpha_sonsuz * c / (4 * b)
mu_hesaplama_2 = mu * (np.array(alpha_hucum) - alpha_zero_lift)
mu_hesaplama_2_rad = np.deg2rad(mu_hesaplama_2)

# FORMÜL KISMI
teta_grafik = []
for i in range(1, 2 * N + 1, 2):
    teta_grafik_iter = np.sin(i * theta) * (1 + (mu * i) / (np.sin(list(theta))))
    teta_grafik.append(teta_grafik_iter)

denklem = np.asarray(teta_grafik)
denklem_transpoz = np.transpose(denklem)
denklem_is_active = np.linalg.inv(denklem_transpoz)

new_matrix = np.matmul(denklem_is_active, mu_hesaplama_2_rad)
alan_deneme_sayilari = np.divide((4 * b), c)

iter0 = (np.sin((1) * theta)) * new_matrix[0] * alan_deneme_sayilari
iter1 = (np.sin((3) * theta)) * new_matrix[1] * alan_deneme_sayilari
iter2 = (np.sin((5) * theta)) * new_matrix[2] * alan_deneme_sayilari
iter3 = (np.sin((7) * theta)) * new_matrix[3] * alan_deneme_sayilari
iter4 = (np.sin((9) * theta)) * new_matrix[4] * alan_deneme_sayilari
iter5 = (np.sin((11) * theta)) * new_matrix[5] * alan_deneme_sayilari
iter6 = (np.sin((13) * theta)) * new_matrix[6] * alan_deneme_sayilari
iter7 = (np.sin((15) * theta)) * new_matrix[7] * alan_deneme_sayilari
iter8 = (np.sin((17) * theta)) * new_matrix[8] * alan_deneme_sayilari
iter9 = (np.sin((19) * theta)) * new_matrix[9] * alan_deneme_sayilari
iter10 = (np.sin((21) * theta)) * new_matrix[10] * alan_deneme_sayilari

# TAŞIMA KATSAYISI HESABI VE DEĞERLERİN OKUNABİLMESİ İÇİN MATRİSE ATANMASI
CL_grafik = iter0 + iter1 + iter2 + iter3 + iter4 + iter5 + iter6 + iter7 + iter8 + iter9 + iter10 # Grafik 'e bastılacak değerler
CL_hesap = (math.pi * AR * new_matrix[0]) # Gerçek hesaplanan CL değeri

rho = (2 * (new_matrix[1] / new_matrix[0])) + (3 * (new_matrix[2] / new_matrix[1])) + (4 * (new_matrix[3] / new_matrix[2])) + (5 * (new_matrix[4] / new_matrix[3])) + (6 * (new_matrix[5] / new_matrix[4])) + (7 * (new_matrix[6] / new_matrix[5]) + (8 * (new_matrix[7] / new_matrix[6]))) + (9 * (new_matrix[8] / new_matrix[7])) + (10 * (new_matrix[9] / new_matrix[8])) + (11 * (new_matrix[10] / new_matrix[9]))

CL_induklenmis_grafik = CL_grafik ** 2 / (math.pi * AR) * rho # Grafik 'e bastılacak değerler
CL_induklenmis_hesap = CL_hesap ** 2 / (math.pi * AR) * rho # Gerçek hesaplanan CL_indüklenmiş değeri

# GRAFİK OLUŞTURMA ALANI
CL1_dizi = np.append(0, CL_grafik)
CL_induklenmis_dizi = np.append(0, CL_induklenmis_grafik)
y_eksen = [b / 2, k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8], k[9], k[10]]

# HESAPLANAN DEĞERLERİ BASTIRMA ALANI
print("Tasima Katsayisi: ", CL_hesap)
print("Induklenmis Surukleme Katsayisi: ", CL_induklenmis_hesap)

# GRAFİK BASTIRMA ALANI
fig = plt.figure("AERODİNAMİK ÖDEVİ [12.GRUP] CL")
plt.plot(y_eksen, CL1_dizi, marker="o")
plt.title("NACA 65(2)-215")
plt.xlabel("Kanat Yarı-Açıklık Konumu (m)")
plt.ylabel("Taşıma Katsayısı")
plt.grid()
plt.show()

fig2 = plt.figure("AERODİNAMİK ÖDEVİ [12.GRUP] CL_induklenmis")
plt.plot(y_eksen, CL_induklenmis_dizi, marker="o")
plt.title("NACA 65(2)-215")
plt.xlabel("Kanat Yarı-Açıklık Konumu (m)")
plt.ylabel("İndüklenmiş Sürükleme Katsayısı")
plt.grid()
plt.show()
