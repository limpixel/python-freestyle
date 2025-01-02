import math
import random
import sys

def mock_input(prompt=""):
    # Fungsi untuk membaca input dari file
    try:
        return input_stream.pop(0)
    except IndexError:
        print("Error: Tidak ada input yang cukup dalam file.")
        sys.exit(1)

# Fungsi dan turunannya
def f_a(x):
    return math.sin(x) - 0.3 * math.exp(x)

def f_b(x):
    return 0.1 * x**3 - 5 * x**2 - x + 4 + math.exp(-x)

def df_a(x):
    return math.cos(x) - 0.3 * math.exp(x)

def df_b(x):
    return 0.3 * x**2 - 10 * x - 1 - math.exp(-x)

# Fungsi dan turunan untuk sistem persamaan nonlinier (persamaan 3)
def f1_system(x, y, h):
    return x - (1 + h**2 * (math.exp(y) * math.sqrt(x) + 3 * x**2))

def f2_system(x, y, h):
    return y - (0.5 + h**2 * math.tan(math.exp(x) + y**2))

def df1_dx(x, y, h):
    return 1 - h**2 * (0.5 * math.exp(y) / math.sqrt(x) + 6 * x)

def df1_dy(x, y, h):
    return -h**2 * math.exp(y) * math.sqrt(x)

def df2_dx(x, y, h):
    return -h**2 * (1 / (math.cos(math.exp(x) + y**2)**2)) * math.exp(x)

def df2_dy(x, y, h):
    return 1 - h**2 * (2 * y / (math.cos(math.exp(x) + y**2)**2))

# Persamaan Van der Waals
def f_vdw(v, P, T, R, a, b):
    return (P + a / v**2) * (v - b) - R * T

def df_vdw(v, P, R, a, b):
    return -2 * a / v**3 * (v - b) + (P + a / v**2)

# Metode Bagidua
def bisection_method(f, a, b, tolerance=1e-6, max_iter=100):
    if a > b:
        a, b = b, a

    iteration = 0
    print(f"{'Iterasi':<10}{'Akar (c)':<15}{'Lebar Selang (b-a)':<20}")
    print("-" * 45)

    while (a - b) < tolerance and iteration < max_iter:
        c = (a + b) / 2
        print(f"{iteration:<10}{c:<15.8f}{(b - a):<20.8f}")
        if f(c) == 0 or (b - a) < tolerance:
            return c, iteration
        iteration += 1
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2, iteration

# Metode Regula Falsi
def regula_falsi_method(f, a, b, tolerance=1e-6, max_iter=100):
    if a > b:
        a, b = b, a
    
    iteration = 0
    print(f"{'Iterasi':<10}{'Akar (c)':<15}{'Lebar Selang (b-a)':<20}{'f(c)':<30}")
    print("-" * 45)

    

    while abs(a - b) > tolerance and iteration < max_iter:
        c = b - (f(b) * (b - a) / (f(b) - f(a)))
        print(f"{iteration:<10}{c:<15.8f}{(b - a):<20.8f}{(f(c)):<30.8f}")
        
        if abs(f(c)) <= tolerance:
            return c, iteration

        iteration += 1
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    
    return c, iteration

# Metode Newton-Raphson
def newton_raphson_method(f, df, x0, tolerance=1e-6, max_iter=100):
    x = x0
    iteration = 0
    print(f"{'Iterasi':<10}{'Xr':<15}{'|X{r+1} - Xr|':<20}")
    print("-" * 45)

    while iteration < max_iter:
        fx = f(x)
        dfx = df(x)

        if dfx == 0:
            raise ValueError("Turunan nol, metode gagal.")

        x_next = x - fx / dfx
        delta = abs(x_next - x)

        print(f"{iteration:<10}{x:<15.8f}{delta:<20.8f}")
        
        if delta < tolerance:
            return x_next, iteration
        
        x = x_next
        iteration += 1
    
    raise ValueError("Metode tidak konvergen.")

# Metode Secant
def secant_method(f, x0, x1, tolerance=1e-6, max_iter=100):
    iteration = 0
    print(f"{'Iterasi':<10}{'Xr':<15}{'|X{r+1} - Xr|':<20}")
    print("-" * 45)

    while iteration < max_iter:
        fx0 = f(x0)
        fx1 = f(x1)

        x_new = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        delta = abs(x_new - x1)    
        
        print(f"{iteration:<10}{x1:<15.8f}{delta:<20.8f}")

        if delta < tolerance:
            return x_new, iteration

        x0, x1 = x1, x_new
        iteration += 1
    
    raise ValueError("Metode tidak konvergen.")

# Metode Newton-Raphson untuk sistem persamaan nonlinier
def newton_raphson_system(h, x0, y0, tolerance, max_iter=100):
    x, y = x0, y0
    print(f"{'Iterasi':<10}{'x':<15}{'y':<15}{'f1(x,y)':<15}{'f2(x,y)':<15}{'delta_x':<15}{'delta_y':<15}")
    for i in range(max_iter):
        f1_val = f1_system(x, y, h)
        f2_val = f2_system(x, y, h)

        df1_dx_val = df1_dx(x, y, h)
        df1_dy_val = df1_dy(x, y, h)
        df2_dx_val = df2_dx(x, y, h)
        df2_dy_val = df2_dy(x, y, h)

        jacobian = [[df1_dx_val, df1_dy_val],
                    [df2_dx_val, df2_dy_val]]

        det = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0]
        if det == 0:
            raise ValueError("Jacobian singular, metode gagal.")

        inv_jacobian = [[jacobian[1][1] / det, -jacobian[0][1] / det],
                        [-jacobian[1][0] / det, jacobian[0][0] / det]]

        delta_x = -(inv_jacobian[0][0] * f1_val + inv_jacobian[0][1] * f2_val)
        delta_y = -(inv_jacobian[1][0] * f1_val + inv_jacobian[1][1] * f2_val)

        x_new = x + delta_x
        y_new = y + delta_y

        if abs(x_new - x) < tolerance and abs(y_new - y) < tolerance:
            return x_new, y_new, i + 1

        x, y = x_new, y_new

    raise ValueError("Metode tidak konvergen dalam iterasi maksimum.")

# Iterasi Gauss-Seidel
def gauss_seidel_system(h, x0, y0, tolerance, max_iter=100):
    x, y = x0, y0
    print(f"{'Iterasi':<10}{'x':<15}{'y':<15}{'delta_x':<15}{'delta_y':<15}")
    for i in range(max_iter):
        x_new = 1 + h**2 * (math.exp(y) * math.sqrt(x) + 3 * x**2)
        y_new = 0.5 + h**2 * math.tan(math.exp(x_new) + y**2)

        if abs(x_new - x) < tolerance and abs(y_new - y) < tolerance:
            return x_new, y_new, i + 1

        x, y = x_new, y_new

    raise ValueError("Metode tidak konvergen dalam iterasi maksimum.")


# Newton-Raphson untuk Volume Gas
def newton_raphson_vdw(P, T, R, a, b, v0, tolerance=1e-6, max_iter=100, show_iterations=False):
    v = v0
    for i in range(max_iter):
        fv = f_vdw(v, P, T, R, a, b)
        dfv = df_vdw(v, P, R, a, b)

        if abs(fv) < tolerance:
            return v, i + 1

        if dfv == 0:
            raise ValueError("Turunan nol, metode gagal.")

        v_new = v - fv / dfv

        if show_iterations:
            print(f"Iterasi {i}: v = {v:.6f}, f(v) = {fv:.6e}")

        if abs(v_new - v) < tolerance:
            return v_new, i + 1

        v = v_new

    raise ValueError("Metode tidak konvergen dalam iterasi maksimum.")

# Fungsi untuk menghitung persamaan fisika
# Fungsi untuk menghitung volume berdasarkan persamaan (a)
def calculate_volume(Vc, d, r, l, theta):
    return Vc +math.pi * d ** 2 * r*(1 - math.cos(theta) + l * (1 - (r/l * math.sin(theta))**2))

# Fungsi untuk persamaan (b) dengan suhu T sebagai variabel
def equation_b(T, Ti, Vi, V, A, B, C, D):
    # Debugging
    print(f"T = {T}, Ti = {Ti}, V = {V}, Vi = {Vi}")
    
    # Validasi input
    if T <= 0 or Ti <= 0:
        raise ValueError("T dan Ti harus positif untuk perhitungan logaritma.")
    if V <= 0 or Vi <= 0:
        raise ValueError("V dan Vi harus positif untuk perhitungan logaritma.")
    
    term1 = A * math.log(T / Ti)
    term2 = B * (T - Ti)
    term3 = 0.5 * C * (T**2 - Ti**2)
    term4 = D * math.log(V / Vi)
    return A * (T - Ti) + B * (T**2 - Ti**2) + C * (T**3 - Ti**3) - D * (V - Vi)

# Fungsi untuk menghitung tekanan berdasarkan persamaan (c)
def calculate_pressure(Pi, Ti, T, Vi, V):
    return Pi * (Vi / V) * (T / Ti)

# Persamaan Gas Ideal
def ideal_gas(P, T, R):
    return R * T / P

# Metode Regula Falsi untuk mencari akar T dari persamaan (b)
def regula_falsi_physic(f, T_lower, T_upper, tolerance=1e-6, max_iter=100, show_iterations=False):
    if show_iterations:
        print(f"{'r':<5}{'a':<12}{'c':<12}{'b':<12}{'f(a)':<12}{'f(c)':<12}{'f(b)':<12}{'Selang Baru':<15}{'Lebarnya':<12}")
        print("-" * 80)

    a, b = T_lower, T_upper
    fa, fb = f(a), f(b)

    for r in range(max_iter):
        c = b - fb * (b - a) / (fb - fa)
        fc = f(c)

        if show_iterations:
            selang_baru = f"[{a:.6f}, {c:.6f}]" if fa * fc < 0 else f"[{c:.6f}, {b:.6f}]"
            lebar = abs(b - a)
            print(f"{r:<5}{a:<12.6f}{c:<12.6f}{b:<12.6f}{fa:<12.6f}{fc:<12.6f}{fb:<12.6f}{selang_baru:<15}{lebar:<12.6f}")

        if abs(fc) < tolerance or abs(b - a) < tolerance:
            return c, r

        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc

    raise ValueError("Metode tidak konvergen.")


# Generate nilai acak
def generate_random_interval(f):
    while True:
        a = random.uniform(-10, 10)
        b = random.uniform(-10, 10)
        if a < b and f(a) * f(b) < 0:
            return a, b

def generate_random_guess():
    return random.uniform(-10, 10)

def right_input(prompt):
    while True:
        try:
            value = float(input(prompt))
            if value == 0:
                print("Nilai tidak boleh 0. Silakan masukkan nilai lain.")
            else:
                return value
        except ValueError:
            print("Input tidak valid. Silakan masukkan angka.")

# Program Utama
def main():
    global input_stream
    input_stream = []
    
    # Cek jika ada file input yang diberikan melalui argument command line
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        try:
            with open(input_file, 'r') as f:
                input_stream = f.read().splitlines()
        except FileNotFoundError:
            print(f"Error: File {input_file} tidak ditemukan.")
            sys.exit(1)
    
    # Gunakan mock_input menggantikan input()
    input_original = __builtins__.input
    __builtins__.input = mock_input
    
    # Jalankan program utama
    try:
        while True:
            print("\n==========================================")
            print("Pilih fungsi:")
            print("1. f(x) = sin(x) - 0.3e^x")
            print("2. f(x) = 0.1x^3 - 5x^2 - x + 4 + e^-x")
            print("3. Hitung 2 persamaan x dan y dengan sistem nonlinear")
            print("4. Hitung persamaan Kimia")
            print("5. Hitung persamaan Fisika")
            print("6. Keluar")
            func_choice = int(input("Masukkan pilihan: "))

            if func_choice == 1:
                f = f_a
                df = df_a
                Four_nonlinear_function(f, df)
            elif func_choice == 2:
                f = f_b
                df = df_b
                Four_nonlinear_function(f, df)
            elif func_choice == 3:
                h = right_input("\nMasukkan nilai h: ")
                x0 = right_input("Masukkan tebakan awal x: ")
                y0 = right_input("Masukkan tebakan awal y: ")
                epsilon = 0.1 * h**4
                nonlinear_function(h, x0, y0, epsilon)
            elif func_choice == 4:
                chemistry_method()
            elif func_choice == 5:
                Pi = right_input("Masukkan tekanan awal Pi (bar): ")
                Ti = right_input("Masukkan temperatur awal Ti (°R): ")
                Vc = right_input("Masukkan volume clearance Vc (m^3): ")
                d = right_input("Masukkan diameter silinder d (m): ")
                r = right_input("Masukkan panjang jari-jari r (m): ")
                l = right_input("Masukkan panjang penghubung l (m): ")
                theta = right_input("Masukkan sudut theta (°): ")
                Vi = right_input("Masukkan volume awal Vi (m^3): ")
                physical_function(Pi, Ti, Vc, d, r, l, theta, Vi)
            elif func_choice == 6:
                print("Terima kasih telah menggunakan program ini.")
                break
            else:
                print("Pilihan tidak valid.")
                return
    finally:
        # Kembalikan input ke fungsi asli
        __builtins__.input = input_original

def Four_nonlinear_function(f, df):
    print("\nPilih metode:")
    print("1. Metode Bagidua (Bisection)")
    print("2. Metode Regula-Falsi")
    print("3. Metode Newton-Raphson")
    print("4. Metode Secant")
    method_choice = int(input("Masukkan pilihan (1, 2, 3, atau 4): "))

    try:
        if method_choice in [1, 2]:  # Metode Bagidua atau Regula Falsi
            a = input("\nMasukkan nilai a (selang awal, tekan Enter untuk random): ")
            b = input("Masukkan nilai b (selang akhir, tekan Enter untuk random): ")
            if not a or not b:  # Jika input kosong
                a, b = generate_random_interval(f)
                print(f"\nNilai acak yang dihasilkan: a = {a}, b = {b}")
            else:
                a, b = float(a), float(b)
            print("\nHasil Iterasi:")
            if method_choice == 1:
                root, iters = bisection_method(f, a, b)
            elif method_choice == 2:
                root, iters = regula_falsi_method(f, a, b)
            print(f"\nAkar ditemukan: {root} setelah {iters} iterasi")

        elif method_choice == 3:  # Metode Newton-Raphson
            x0 = input("\nMasukkan nilai tebakan awal (x0, tekan Enter untuk random): ")
            if not x0:
                x0 = generate_random_guess()
                print(f"\nNilai acak yang dihasilkan: x0 = {x0}")
            else:
                x0 = float(x0)
            print("\nHasil Iterasi:")
            root, iters = newton_raphson_method(f, df, x0)
            print(f"\nAkar ditemukan: {root} setelah {iters} iterasi")

        elif method_choice == 4:  # Metode Secant
            x0 = input("\nMasukkan nilai tebakan awal pertama (x0, tekan Enter untuk random): ")
            x1 = input("Masukkan nilai tebakan awal kedua (x1, tekan Enter untuk random): ")
            if not x0 or not x1:
                x0, x1 = generate_random_guess(), generate_random_guess()
                print(f"\nNilai acak yang dihasilkan: x0 = {x0}, x1 = {x1}")
            else:
                x0, x1 = float(x0), float(x1)
            print("\nHasil Iterasi:")
            root, iters = secant_method(f, x0, x1)
            print(f"\nAkar ditemukan: {root} setelah {iters} iterasi")

        else:
            print("Pilihan metode tidak valid.")
    except ValueError as e:
        print(f"\nError: {e}")

def nonlinear_function(h, x0, y0, epsilon):    
    print("\nPilih metode:")
    print("1. Metode Newton-Raphson")
    print("2. Iterasi Gauss-Seidel")
    method_choice = int(input("Masukkan pilihan (1 atau 2): "))

    if method_choice == 1:
        x, y, iters = newton_raphson_system(h, x0, y0, epsilon)
        print(f"\nHasil: x = {x}, y = {y}, iterasi = {iters}")
    elif method_choice == 2:
        x, y, iters = gauss_seidel_system(h, x0, y0, epsilon)
        print(f"\nHasil: x = {x}, y = {y}, iterasi = {iters}")
    else:
        print("Pilihan metode tidak valid.")

def chemistry_method():
    P = 70  # bar
    T = 215  # Kelvin
    R = 0.08314  # bar m3/(kg mol K)
    a = 1.463  # bar m6/(kg mol)^2
    b = 0.0394  # m3/kg

    # Tebakan awal untuk Newton-Raphson
    v0 = 0.1  # m3/kg (nilai awal mendekati)

    try:
        # Hitung dengan Van der Waals
        v_vdw, iter_vdw = newton_raphson_vdw(P, T, R, a, b, v0, show_iterations=True)
        print(f"Volume (Van der Waals): {v_vdw:.6f} m3/kg setelah {iter_vdw} iterasi")

        # Hitung dengan gas ideal
        v_ideal = ideal_gas(P, T, R)
        print(f"Volume (Gas Ideal): {v_ideal:.6f} m3/kg")

        # Perbandingan
        error = abs(v_vdw - v_ideal) / v_ideal * 100
        print(f"Perbedaan relatif: {error:.2f}%")
    except ValueError as e:
        print(f"Error: {e}")

def physical_function(Pi, Ti, Vc, d, r, l, theta,Vi):
    # Konstanta untuk persamaan (b)
    A = 0.15787
    B = 0.51001e-4
    C = 0.74171e-8
    D = 0.6855

    # Hitung volume berdasarkan persamaan (a)
    V = calculate_volume(Vc, d, r, l, theta)
    print(f"Volume V berdasarkan persamaan (a): {V:.6f} m^3")

    # Cari temperatur T menggunakan metode Regula Falsi pada persamaan (b)
    try:
        T_lower = Ti  # Tebakan awal bawah (mendekati Ti)
        T_upper = Ti + 100  # Tebakan awal atas
        f = lambda T: equation_b(T, Ti, Vi, V, A, B, C, D)
        T, iter_count = regula_falsi_physic(f, T_lower, T_upper, show_iterations=True)
        print(f"Temperatur T berdasarkan persamaan (b): {T:.6f} °R setelah {iter_count} iterasi")
    except ValueError as e:
        print(f"Error dalam menghitung T: {e}")
        return
    
    # Hitung tekanan berdasarkan persamaan (c)
    P = calculate_pressure(Pi, Ti, T, Vi, V)
    print(f"Tekanan P berdasarkan persamaan (c): {P:.6f} bar")

if __name__ == "__main__":
    main()