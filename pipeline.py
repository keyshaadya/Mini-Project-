import csv
import os
import matplotlib.pyplot as plt

# ==========================================================
# 1. MEMBACA FILE FASTA & DISIMPAN DALAM LIST
# ==========================================================
def baca_fasta(nama_file):
    # Validasi apakah file sekuens.fasta berada di folder yang sama
    if not os.path.exists(nama_file):
        print(f"[EROR] File '{nama_file}' tidak ditemukan di folder aktif!")
        print("Pastikan file FASTA hasil download dari NCBI sudah dipindahkan")
        print(f"ke folder proyek dan diberi nama '{nama_file}'")
        return None
        
    sekuens_list = []  # Struktur Data List sesuai instruksi soal
    current_id = ""
    current_seq = ""
    
    with open(nama_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    sekuens_list.append({"id": current_id, "sekuens": current_seq})
                current_id = line[1:]  # Mengambil ID tanpa tanda '>'
                current_seq = ""       
            else:
                current_seq += line
                
        if current_id:
            sekuens_list.append({"id": current_id, "sekuens": current_seq})
            
    return sekuens_list

# ===================================================================
# 2. MENGHITUNG FREKUENSI NUKLEOTIDA & GC CONTENT (DICTIONARY)
# ===================================================================
def analisis_pipeline(sekuens_list):
    hasil_analisis = [] 
    
    for item in sekuens_list:
        seq = item["sekuens"].upper() 
        
        # Menggunakan Dictionary sesuai instruksi soal untuk menghitung frekuensi
        frekuensi = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for nukleotida in seq:
            if nukleotida in frekuensi:
                frekuensi[nukleotida] += 1
                
        total_basa = len(seq)
        gc_count = frekuensi['G'] + frekuensi['C']
        gc_content = (gc_count / total_basa) * 100 if total_basa > 0 else 0
        
        hasil_analisis.append({
            "id": item["id"],
            "A": frekuensi['A'],
            "T": frekuensi['T'],
            "G": frekuensi['G'],
            "C": frekuensi['C'],
            "GC_Content": round(gc_content, 2)
        })
        
    return hasil_analisis

# ===================================================================
# 3. EKSEKUSI UTAMA PIPELINE
# ===================================================================
NAMA_FILE_FASTA = "sekuens.fasta" 

# 1. Jalankan pembacaan file FASTA
data_mentah = baca_fasta(NAMA_FILE_FASTA)

if data_mentah:
    # 2. Jalankan analisis frekuensi & GC Content
    hasil_proses = analisis_pipeline(data_mentah)

    # 4. MENGURUTKAN BERDASARKAN GC CONTENT (SORTING)
    hasil_terurut = sorted(hasil_proses, key=lambda x: x['GC_Content'], reverse=True)

    # 5. MENAMPILKAN 3 SEKUENS TERBAIK DI TERMINAL VS CODE
    print("\n=== 3 SEKUENS TERBAIK BERDASARKAN GC CONTENT ===")
    tiga_terbaik = hasil_terurut[:3] # Slicing list untuk mengambil 3 data teratas
    for i, data in enumerate(tiga_terbaik, 1):
        print(f"{i}. ID: {data['id'][:50]}... \n   GC Content: {data['GC_Content']}%\n")

    # 6. MENULISKAN HASIL KE FILE CSV
    with open("hasil_analisis.csv", "w", newline="") as csv_file:
        fieldnames = ["id", "A", "T", "G", "C", "GC_Content"]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(hasil_proses)
    print("[INFO] Sukses! File 'hasil_analisis.csv' berhasil dibuat di foldermu.")

    # ===================================================================
    # 7. VISUALISASI GRAFIK HASIL BERDASARKAN NILAI GC (MATPLOTLIB)
    # ===================================================================
    # Mengambil ringkasan ID agar tampilan grafik tidak bertumpuk teksnya
    ids_singkat = [data['id'].split()[0][:15] for data in hasil_proses]
    gc_values = [data['GC_Content'] for data in hasil_proses]

    plt.figure(figsize=(9, 5))
    # Membuat grafik batang dengan warna yang kontras dan elegan
    colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c'][:len(gc_values)]
    bars = plt.bar(ids_singkat, gc_values, color=colors, edgecolor='black', width=0.5)

    plt.ylabel('GC Content (%)', fontsize=11, fontweight='bold')
    plt.xlabel('Sekuens ID (Ringkas)', fontsize=11, fontweight='bold')
    plt.title('Visualisasi Analisis GC Content - Pipeline Gen PETase', fontsize=13, fontweight='bold', pad=15)
    plt.ylim(0, 100)
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    # Menambahkan nilai persentase tepat di atas setiap batang grafik
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 2,
                 f'{height}%',
                 ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()
    
    # Otomatis menyimpan grafik menjadi file gambar PNG di foldermu
    plt.savefig('grafik_gc_content.png', dpi=300) 
    print("[INFO] Sukses! Gambar 'grafik_gc_content.png' berhasil disimpan.")
    
    # Menampilkan jendela grafik popup di laptopmu
    plt.show()