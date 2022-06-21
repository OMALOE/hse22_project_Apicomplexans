# Apicomplexans
<h3> Участники </h3>

| Участник  | Род | Git | 
|---|---|---|
| Архипов Дмитрий  |Eimeria (mitis, praecox, maxima, acervulina, necatrix) | https://github.com/OMALOE/hse22_project_bioinf |
| Волгунов Фёдор   | Theileria, Babesia |https://github.com/hq43et28ms9z/hse22_project |
| Дыбовский Леонид  | Plasmodium | https://github.com/LeonidDybovskij/hse_project_plasmodium |
| Марусев Егор  | Cryptosporidium    | https://github.com/egormarusev/hse22_project_bio.git |
| Галкина Мария  | Babesia| https://github.com/galkinamariia/hse22_project |
| Пирожкина Мария  |Porospora, Gregarina niphandrodes, Nephromyces, Haemoproteus  | https://github.com/Pirozhok1967/hse22_project |
| Тарасова Мария  |Plasmodium( relictum, reichenowi, berghei ANKA), Besnoitia, Neospora   | https://github.com/MariaTar7/hse22_project/ |
| Шмелев Алексей |Plasmodium (chabaudi, fragile, gallinaceum, coatneyi, malariae)| https://github.com/alexeyshmelev/hse22_project |
| Жукова Юлия  | Eimeria   (brunetti,  tenella, maxima, acervulina, necatrix)    | https://github.com/ZhukovaJul/hse22_project.git |
| Смирнов Георгий  |Toxoplasma, Plasmodium(sp. gorilla clade G2, cynomolgi strain B, vinckei, inui San Antonio 1) | пока нет |

<h3> Взятые виды: </h3>

##### *Plasmodium*
 * *falciparum*
 * *vivax*
 * *gaboni*
 * *knowlesi*
 * *yoelii*
 * *relictum*
 * *reichenowi*
 * *berghei ANKA*
 * *chabaudi*
 * *fragile*
 * *gallinaceum*
 * *coatneyi*
 * *malariae*
 * *vinckei*
 * *sp. gorilla clade G2*
 * *inui San Antonio 1*
 * *cynomolgi strain B*
##### *Eimeria*
 * *mitis*
 * *necatrix*
 * *praecox*
 * *maxima*
 * *acervulina*
 * *brunetti*
 * *tenella*
##### *Theileria*
 * *annulata*
 * *equi strain WA*
 * *orientalis strain Shintoku*
 * *parva strain Muguga*
##### *Babesia*
 * *bovis T2Bo*
 * *microti strain RI*
 * *bigemina*
 * *sp. Xinjiang*
 * *ovata*
##### *Cryptosporidium*
 * *parvum Iowa II*
 * *hominis TU502*
 * *muris RN66*
 * *ubiquitum*
 * *felis*
##### *Gregarina*
 * *niphandrodes*
##### *Porospora*
 * *cf. gigantea A*
 * *cf. gigantea B*
##### *Besnoitia*
 * *besnoiti*
##### *Neospora*
 * *caninum Liverpool*
##### *Toxoplasma*
 * *gondii ME49*


<h3> Создание ортологичных кластеров внутри таксона </h3> 
Для создания кластеров был запущен proteinortho на 44 геномах, скачанных с помощью скрипта на bash (папка download_files).
Результат - в файле full_gca_44.tsv

<h3> Статистика по кластерам  </h3> 
Для дополнения таблицы статистикой был использован скрипт на python (файл group_clusters.py, папка src)
Всего было построено 26462 кластера. Их распределение представлено на двух картинках ниже.

Распределение числа геномов, содержащих белки в кластере (количество ортологов):

![gcf_gca_44_clust_by_number_of_genoms](https://user-images.githubusercontent.com/60808642/174497142-7afd4c5f-5239-4322-ba56-fbc879cc2ae9.png)

Распределение числа белков на кластер (с учетом числа паралогов внутри каждого вида):

![gcf_gca_44_clust_by_number_of_proteins](https://user-images.githubusercontent.com/60808642/174497029-1004d50f-6a03-4200-b6c4-a88ddbb0a272.png)

Как видно из графиков, большинство кластеров охватывает всего несколько видов. Поэтому для дальнейшей работы была взята отсечка в 22 генома на кластер. Таких кластеров оказалось 2363, их распределение показано на рисунках ниже.

Распеределение числа геномов, содержащих белки в кластере (количество ортологов):

![full_gca_cluster_g_stat_22_plus](https://user-images.githubusercontent.com/60808642/174497118-0f541621-4200-4005-905a-33f4077524ed.png)

Распределение числа белков на кластер (с учетом числа паралогов внутри каждого вида):

![full_gca_cluster_p_stat_22_plus](https://user-images.githubusercontent.com/60808642/174497123-52c7d0ce-a3b0-4509-bb29-ebb24ddcb8b3.png)

Нужно отметить, что среди этих 2363 кластеров был 51 кластер, содержащий белки из всех 44 видов.
Полная статистика по кластерам - в файле full_gca_cluster_stat_with_function.tsv

<h3> Выбор кластеров </h3> 

Дополнение таблицы числами Z-SCORE осуществлялось с помощью файла group.py, работающего со скачанными для 44 видов файлами fna и feature_table.txt (папка src). Результат - в файле clusters_with_z_score.tsv.

<h3> Тепловая карта  </h3> 

Серым обозначены ячейки, имеющие score < 500, а чёрным - ячейки, не имеющие score совсем (т.е. нет пересечения с генами, в которых есть Z-DNA). Также стоит отметить, что на данной картинке изображены первые 20 из 2363 найденных кластеров (сортировка по убыванию была произведена по среднему score в строке).

![clusters](https://user-images.githubusercontent.com/60858323/174839838-b834b3a4-134f-46d0-a298-a8e6a5a0debf.jpg)


<h3> Функции кластеров </h3> 

<h3> Визуализация расположение участков Z-DNA относительно генов </h3> 

<h3> Выравнивания </h3> 
