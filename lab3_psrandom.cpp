/// \file lab3_psrandom.cpp
/// \brief Лабораторная работа №3. Исследование генераторов псевдослучайных чисел.
/// \details
/// Программа реализует три пользовательских генератора псевдослучайных чисел,
/// сравнивает их статистические характеристики, выполняет упрощённые NIST-подобные
/// тесты и измеряет скорость генерации. Результаты сохраняются в CSV-файлы
/// для последующего построения графиков и включения в отчёт Doxygen.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace std;

/// \brief Количество выборок для статистического анализа.
const int NUM_SAMPLES = 20;

/// \brief Размер одной выборки.
const int SAMPLE_SIZE = 1000;

/// \brief Размер диапазона случайных чисел.
const int RANGE_SIZE = 10000;

/// \brief Количество интервалов для критерия хи-квадрат.
const int K = 20;

/// \brief Размер одного интервала для критерия хи-квадрат.
const int BIN_SIZE = RANGE_SIZE / K;

/// \brief Критическое значение критерия хи-квадрат при 19 степенях свободы и alpha = 0.05.
const double CHI2_CRIT = 30.1435;

/// \brief Количество битов для NIST-подобных тестов.
const int BIT_COUNT = 1000000;

/// \brief Циклический сдвиг 32-битного числа влево.
/// \param x Исходное значение.
/// \param r Количество бит для сдвига.
/// \return Результат циклического сдвига.
inline uint32_t rotl32(uint32_t x, unsigned r) {
    return (x << r) | (x >> (32 - r));
}

/// \brief Генератор на основе линейного конгруэнтного метода с перемешиванием Bays-Durham.
/// \details
/// Базовая последовательность получается линейным конгруэнтным генератором.
/// Далее используется таблица перемешивания из 32 элементов для уменьшения корреляций.
class LCG_BaysDurham {
    static const uint64_t A = 6364136223846793005ULL;   ///< Множитель LCG.
    static const uint64_t C = 1442695040888963407ULL;   ///< Приращение LCG.
    static const int table_size = 32;                  ///< Размер таблицы перемешивания.

    uint64_t state;                                     ///< Текущее состояние генератора.
    array<uint64_t, table_size> table{};               ///< Таблица перемешивания.
    uint64_t last;                                      ///< Последнее использованное значение.

    /// \brief Выполняет один шаг базового линейного конгруэнтного генератора.
    /// \return Следующее 64-битное значение.
    uint64_t rawLCG() {
        state = A * state + C;
        return state;
    }

public:
    /// \brief Конструктор генератора.
    /// \param seed Начальное зерно.
    explicit LCG_BaysDurham(uint64_t seed = 123456789ULL) {
        if (seed != 0)
            state = seed;
        else
            state = 123456789ULL;
        last = 0;
        for (int i = 0; i < table_size; ++i) table[i] = rawLCG();
        last = rawLCG();
    }

    /// \brief Переинициализация генератора новым зерном.
    /// \param seedValue Новое зерно.
    void seed(uint64_t seedValue) {
        state = seedValue ? seedValue : 123456789ULL;
        for (int i = 0; i < table_size; ++i) table[i] = rawLCG();
        last = rawLCG();
    }

    /// \brief Генерация 64-битного псевдослучайного числа.
    /// \return Следующее 64-битное значение.
    uint64_t next64() {
        uint64_t x = rawLCG();
        size_t j = static_cast<size_t>((last >> 59) & (table_size - 1));
        uint64_t y = table[j];
        table[j] = x;
        last = y;

        y ^= (y >> 27);
        y *= 0x3C79AC492BA7B653ULL;
        y ^= (y >> 33);
        y *= 0x1C69B3F74AC4AE35ULL;
        y ^= (y >> 27);
        return y;
    }

    /// \brief Генерация 32-битного числа.
    /// \return Следующее 32-битное значение.
    uint32_t next() {
        return static_cast<uint32_t>(next64());
    }
};

/// \brief Генератор xorshift128 с дополнительным перемешиванием.
/// \details
/// Используется схема xorshift128, после которой результат дополнительно
/// подвергается смешиванию для улучшения статистических свойств.
class Xorshift128_MLT {
    uint64_t s[2]; ///< Внутреннее состояние генератора.

public:
    /// \brief Конструктор генератора.
    /// \param seed Начальное зерно.
    explicit Xorshift128_MLT(uint64_t seed = 123456789ULL) {
        this->seed(seed);
    }

    /// \brief Переинициализация генератора.
    /// \param seedValue Новое зерно.
    void seed(uint64_t seedValue) {
        if (seedValue == 0) seedValue = 123456789ULL;
        s[0] = seedValue;
        s[1] = seedValue ^ 0x9e3779b97f4a7c15ULL;
        for (int i = 0; i < 16; ++i) next64();
    }

    /// \brief Генерация 64-битного числа.
    /// \return Следующее 64-битное значение.
    uint64_t next64() {
        uint64_t s1 = s[0];
        uint64_t s0 = s[1];
        s[0] = s0;
        s1 ^= s1 << 23;
        s1 ^= s1 >> 17;
        s1 ^= s0;
        s1 ^= s0 >> 26;
        s[1] = s1;

        uint64_t result = s1 + s0;
        result ^= (result >> 30);
        result *= 0xbf58476d1ce4e5b9ULL;
        result ^= (result >> 27);
        result *= 0x94d049bb133111ebULL;
        result ^= (result >> 31);
        return result;
    }

    /// \brief Генерация 32-битного числа.
    /// \return Следующее 32-битное значение.
    uint32_t next() {
        return static_cast<uint32_t>(next64());
    }
};

/// \brief Аддитивный лаговый рекуррентный генератор с вращением и XOR.
/// \details
/// Генератор использует кольцевой буфер, аддитивную схему с двумя лагами,
/// циклические сдвиги и XOR-перемешивание.
class AdditiveLFG_RotateXor {
    static const int LAG_A = 55;   ///< Основной лаг.
    static const int LAG_B = 24;   ///< Дополнительный лаг.

    array<uint32_t, LAG_A> buffer{}; ///< Кольцевой буфер значений.
    int index;                       ///< Текущая позиция в буфере.

    /// \brief Дополнительное перемешивание 32-битного значения.
    /// \param x Исходное значение.
    /// \return Перемешанное значение.
    static uint32_t scramble(uint32_t x) {
        x = rotl32(x, 7) ^ (x >> 3);
        x += 0x9e3779b9U;
        x ^= rotl32(x, 11);
        x ^= (x >> 16);
        return x;
    }

    /// \brief Перезаполнение буфера новыми значениями.
    void refill() {
        for (int i = 0; i < LAG_A; ++i) {
            int j = (i + LAG_B) % LAG_A;
            uint32_t v = buffer[i] + buffer[j];
            v ^= rotl32(buffer[j], 9);
            v = scramble(v);
            buffer[i] = v;
        }
        index = 0;
    }

public:
    /// \brief Конструктор генератора.
    /// \param seed Начальное зерно.
    explicit AdditiveLFG_RotateXor(uint32_t seed = 123456789U) {
        this->seed(seed);
    }

    /// \brief Переинициализация генератора.
    /// \param seedValue Новое зерно.
    void seed(uint32_t seedValue) {
        uint32_t x = seedValue ? seedValue : 123456789U;
        for (int i = 0; i < LAG_A; ++i) {
            x = 1664525U * x + 1013904223U;
            buffer[i] = scramble(x + static_cast<uint32_t>(i * 2654435761U));
        }
        index = LAG_A;
        for (int k = 0; k < 8; ++k) refill();
    }

    /// \brief Генерация 32-битного значения из буфера.
    /// \return Следующее 32-битное значение.
    uint32_t next32() {
        if (index >= LAG_A) refill();
        return buffer[index++];
    }

    /// \brief Генерация 32-битного значения.
    /// \return Следующее 32-битное значение.
    uint32_t next() {
        return next32();
    }
};

/// \brief Структура статистик для одной выборки.
struct SampleStats {
    double mean;     ///< Среднее значение.
    double stddev;   ///< Стандартное отклонение.
    double cv;       ///< Коэффициент вариации.
    double chi2;     ///< Статистика хи-квадрат.
    bool uniform;    ///< Результат проверки на равномерность.
};

/// \brief Вычисляет статистики для выборки случайных чисел.
/// \param sample Выборка целых чисел.
/// \return Заполненная структура статистик.
SampleStats computeStats(const vector<int>& sample) {
    int n = static_cast<int>(sample.size());
    double sum = 0.0;
    for (int x : sample) sum += x;
    double mean = sum / n;

    double sq = 0.0;
    for (int x : sample) sq += (x - mean) * (x - mean);
    double stddev = sqrt(sq / (n - 1));
    double cv = (mean != 0.0) ? (stddev / mean) : 0.0;

    int bins[K] = { 0 };
    for (int x : sample) {
        int bin = x / BIN_SIZE;
        if (bin >= K) bin = K - 1;
        ++bins[bin];
    }

    double expected = double(n) / K;
    double chi2 = 0.0;
    for (int i = 0; i < K; ++i) {
        double diff = bins[i] - expected;
        chi2 += (diff * diff) / expected;
    }

    bool uniform = (chi2 <= CHI2_CRIT);
    return { mean, stddev, cv, chi2, uniform };
}

/// \brief Генерирует выборку случайных чисел заданного размера.
/// \tparam Gen Тип генератора.
/// \param gen Экземпляр генератора.
/// \param size Размер выборки.
/// \param range Размер диапазона значений.
/// \return Вектор случайных чисел.
template<typename Gen>
vector<int> generateSample(Gen& gen, int size, int range) {
    vector<int> sample;
    sample.reserve(size);
    for (int i = 0; i < size; ++i)
        sample.push_back(static_cast<int>(gen.next() % range));
    return sample;
}

/// \brief Аппроксимация дополнительной функции ошибок.
/// \param x Аргумент.
/// \return Значение erfc(x).
double my_erfc(double x) {
    double t = 1.0 / (1.0 + 0.5 * fabs(x));
    double tau = t * exp(-x * x - 1.26551223 + t * (1.00002368 + t * (0.37409196 +
        t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 +
            t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return (x >= 0) ? tau : 2.0 - tau;
}

/// \brief Вычисляет p-value для статистики хи-квадрат.
/// \param chi2 Значение статистики.
/// \param df Число степеней свободы.
/// \return Приближённое p-value.
double p_value_from_chi2(double chi2, int df) {
    double z = (sqrt(2.0 * chi2) - sqrt(2.0 * df - 1.0));
    return my_erfc(fabs(z) / sqrt(2.0));
}

/// \brief NIST-подобный монобитный тест.
/// \param bits Последовательность битов.
/// \return p-value теста.
double nist_monobit(const vector<uint8_t>& bits) {
    int n = static_cast<int>(bits.size());
    if (n < 100) return -1;
    int sum = 0;
    for (auto b : bits) sum += b ? 1 : -1;
    double s = sum / sqrt(static_cast<double>(n));
    return my_erfc(fabs(s) / sqrt(2.0));
}

/// \brief NIST-подобный тест частоты в блоках.
/// \param bits Последовательность битов.
/// \param M Размер блока.
/// \return p-value теста.
double nist_block_frequency(const vector<uint8_t>& bits, int M = 128) {
    int n = static_cast<int>(bits.size());
    if (n < M) return -1;
    int N = n / M;
    double chi2 = 0.0;
    for (int i = 0; i < N; ++i) {
        int ones = 0;
        for (int j = 0; j < M; ++j)
            if (bits[i * M + j]) ++ones;
        double pi = static_cast<double>(ones) / M;
        chi2 += (pi - 0.5) * (pi - 0.5);
    }
    chi2 *= 4.0 * M;
    return p_value_from_chi2(chi2, N);
}

/// \brief NIST-подобный тест серий.
/// \param bits Последовательность битов.
/// \return p-value теста.
double nist_runs(const vector<uint8_t>& bits) {
    int n = static_cast<int>(bits.size());
    if (n < 100) return -1;
    int ones = 0;
    for (auto b : bits) ones += b;
    double pi = static_cast<double>(ones) / n;
    if (fabs(pi - 0.5) >= 2.0 / sqrt(static_cast<double>(n))) return 0.0;
    int runs = 1;
    for (int i = 1; i < n; ++i)
        if (bits[i] != bits[i - 1]) ++runs;
    double num = runs - 2.0 * n * pi * (1.0 - pi);
    double den = 2.0 * sqrt(2.0 * n) * pi * (1.0 - pi);
    double z = num / den;
    return my_erfc(fabs(z) / sqrt(2.0));
}

/// \brief NIST-подобный тест на длину наибольшей серии единиц.
/// \param bits Последовательность битов.
/// \return p-value теста.
double nist_longest_run(const vector<uint8_t>& bits) {
    int n = static_cast<int>(bits.size());
    if (n < 128) return -1;

    const int M = 10000;
    int N = n / M;
    if (N == 0) return -1;
    const int CAT = 7;
    vector<int> v(CAT, 0);

    for (int i = 0; i < N; ++i) {
        int longest = 0, current = 0;
        for (int j = 0; j < M; ++j) {
            if (bits[i * M + j]) ++current;
            else {
                if (current > longest) longest = current;
                current = 0;
            }
        }
        if (current > longest) longest = current;

        if (longest <= 10) v[0]++;
        else if (longest == 11) v[1]++;
        else if (longest == 12) v[2]++;
        else if (longest == 13) v[3]++;
        else if (longest == 14) v[4]++;
        else if (longest == 15) v[5]++;
        else v[6]++;
    }

    double pi[CAT] = { 0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727 };
    double chi2 = 0.0;
    for (int i = 0; i < CAT; ++i) {
        double expected = N * pi[i];
        chi2 += (v[i] - expected) * (v[i] - expected) / expected;
    }
    return p_value_from_chi2(chi2, CAT - 1);
}

/// \brief NIST-подобный тест кумулятивных сумм.
/// \param bits Последовательность битов.
/// \param forward Направление обхода последовательности.
/// \return p-value теста.
double nist_cusum(const vector<uint8_t>& bits, bool forward = true) {
    int n = static_cast<int>(bits.size());
    if (n < 100) return -1;
    vector<int> X(n, 0);
    int sum = 0;
    if (forward) {
        for (int i = 0; i < n; ++i) {
            sum += bits[i] ? 1 : -1;
            X[i] = sum;
        }
    }
    else {
        for (int i = n - 1; i >= 0; --i) {
            sum += bits[i] ? 1 : -1;
            X[n - 1 - i] = sum;
        }
    }
    int maxZ = 0;
    for (auto z : X) if (abs(z) > maxZ) maxZ = abs(z);
    double z_stat = maxZ / sqrt(static_cast<double>(n));
    return my_erfc(z_stat / sqrt(2.0));
}

/// \brief Обёртка над стандартным генератором Mersenne Twister.
/// \details Используется для сравнения скорости с пользовательскими генераторами.
class StandardMT {
    mt19937_64 gen;                               ///< Внутренний генератор.
    uniform_int_distribution<uint32_t> dist;     ///< Распределение на диапазоне uint32_t.
public:
    /// \brief Конструктор генератора.
    /// \param seed Начальное зерно.
    explicit StandardMT(uint64_t seed = 12345)
        : gen(seed), dist(0, numeric_limits<uint32_t>::max()) {
    }

    /// \brief Переинициализация генератора.
    /// \param seedValue Новое зерно.
    void seed(uint64_t seedValue) { gen.seed(seedValue); }

    /// \brief Генерация 32-битного значения.
    /// \return Следующее псевдослучайное число.
    uint32_t next() { return dist(gen); }
};

/// \brief Генерирует битовую последовательность.
/// \tparam Gen Тип генератора.
/// \param gen Экземпляр генератора.
/// \param numBits Количество бит.
/// \return Вектор битов.
template<typename Gen>
vector<uint8_t> generateBits(Gen& gen, int numBits) {
    vector<uint8_t> bits;
    bits.reserve(numBits);
    for (int i = 0; i < numBits; ++i) {
        uint32_t r = gen.next();
        bits.push_back((r >> 31) & 1U);
    }
    return bits;
}

/// \brief Измеряет время генерации заданного количества чисел.
/// \tparam Gen Тип генератора.
/// \param gen Экземпляр генератора.
/// \param count Количество генерируемых значений.
/// \return Время генерации в миллисекундах.
template<typename Gen>
double measureGenerationTime(Gen& gen, int count) {
    auto t0 = chrono::high_resolution_clock::now();
    for (int i = 0; i < count; ++i) {
        volatile uint32_t sink = gen.next();
        (void)sink;
    }
    auto t1 = chrono::high_resolution_clock::now();
    return chrono::duration<double, milli>(t1 - t0).count();
}

/// \brief Точка входа в программу.
/// \details
/// Программа:
/// - вычисляет описательные статистики для 20 выборок каждого генератора;
/// - сохраняет их в файл output/random_stats.csv;
/// - выполняет NIST-подобные тесты и сохраняет результаты в output/nist_results.csv;
/// - измеряет скорость генерации и сохраняет её в output/speed_results.csv.
int main() {
    auto now = chrono::high_resolution_clock::now().time_since_epoch().count();
    uint32_t baseSeed = static_cast<uint32_t>(now);

    FILE* f = fopen("output/random_stats.csv", "w");
    fprintf(f, "generator,sample,mean,stddev,cv,chi2,uniform\n");

    printf("LCG_BaysDurham results:\n");
    for (int s = 0; s < NUM_SAMPLES; ++s) {
        LCG_BaysDurham gen(baseSeed + s * 1000);
        auto sample = generateSample(gen, SAMPLE_SIZE, RANGE_SIZE);
        auto stats = computeStats(sample);
        printf("Sample %2d: mean=%7.1f  stddev=%7.1f  cv=%6.4f  chi2=%7.3f  %s\n",
            s + 1, stats.mean, stats.stddev, stats.cv, stats.chi2,
            stats.uniform ? "PASS" : "FAIL");
        fprintf(f, "LCG_BaysDurham,%d,%.1f,%.1f,%.4f,%.3f,%d\n",
            s + 1, stats.mean, stats.stddev, stats.cv, stats.chi2, stats.uniform);
    }

    printf("\nXorshift128_MLT results:\n");
    for (int s = 0; s < NUM_SAMPLES; ++s) {
        Xorshift128_MLT gen(baseSeed + s * 2000);
        auto sample = generateSample(gen, SAMPLE_SIZE, RANGE_SIZE);
        auto stats = computeStats(sample);
        printf("Sample %2d: mean=%7.1f  stddev=%7.1f  cv=%6.4f  chi2=%7.3f  %s\n",
            s + 1, stats.mean, stats.stddev, stats.cv, stats.chi2,
            stats.uniform ? "PASS" : "FAIL");
        fprintf(f, "Xorshift128_MLT,%d,%.1f,%.1f,%.4f,%.3f,%d\n",
            s + 1, stats.mean, stats.stddev, stats.cv, stats.chi2, stats.uniform);
    }

    printf("\nAdditiveLFG_RotateXor results:\n");
    for (int s = 0; s < NUM_SAMPLES; ++s) {
        AdditiveLFG_RotateXor gen(baseSeed + s * 3000);
        auto sample = generateSample(gen, SAMPLE_SIZE, RANGE_SIZE);
        auto stats = computeStats(sample);
        printf("Sample %2d: mean=%7.1f  stddev=%7.1f  cv=%6.4f  chi2=%7.3f  %s\n",
            s + 1, stats.mean, stats.stddev, stats.cv, stats.chi2,
            stats.uniform ? "PASS" : "FAIL");
        fprintf(f, "AdditiveLFG_RotateXor,%d,%.1f,%.1f,%.4f,%.3f,%d\n",
            s + 1, stats.mean, stats.stddev, stats.cv, stats.chi2, stats.uniform);
    }
    fclose(f);
    printf("\nDescriptive statistics saved to output/random_stats.csv\n");

    FILE* nist_f = fopen("output/nist_results.csv", "w");
    fprintf(nist_f, "generator,monobit,block_freq,runs,longest_run,cusum_forward,cusum_backward\n");

    vector<string> generators = { "LCG_BaysDurham", "Xorshift128_MLT", "AdditiveLFG_RotateXor" };
    for (int g = 0; g < 3; ++g) {
        vector<uint8_t> bits;
        uint32_t seed = baseSeed + g * 5000;
        if (g == 0) {
            LCG_BaysDurham gen(seed);
            bits = generateBits(gen, BIT_COUNT);
        }
        else if (g == 1) {
            Xorshift128_MLT gen(seed);
            bits = generateBits(gen, BIT_COUNT);
        }
        else {
            AdditiveLFG_RotateXor gen(seed);
            bits = generateBits(gen, BIT_COUNT);
        }

        double mono = nist_monobit(bits);
        double block = nist_block_frequency(bits, 128);
        double runs = nist_runs(bits);
        double longest = nist_longest_run(bits);
        double cusum_f = nist_cusum(bits, true);
        double cusum_b = nist_cusum(bits, false);

        printf("\n%s NIST-like results (p-values):\n", generators[g].c_str());
        printf("Monobit:         %.6f\n", mono);
        printf("Block Frequency: %.6f\n", block);
        printf("Runs:            %.6f\n", runs);
        printf("Longest Run:     %.6f\n", longest);
        printf("CUSUM forward:   %.6f\n", cusum_f);
        printf("CUSUM backward:  %.6f\n", cusum_b);

        fprintf(nist_f, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            generators[g].c_str(), mono, block, runs, longest, cusum_f, cusum_b);
    }
    fclose(nist_f);
    printf("\nNIST-style results saved to output/nist_results.csv\n");

    const vector<int> SIZES = { 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000 };
    FILE* speed_f = fopen("output/speed_results.csv", "w");
    fprintf(speed_f, "size,LCG_BaysDurham_ms,Xorshift128_MLT_ms,AdditiveLFG_RotateXor_ms,StandardMT_ms\n");

    printf("\n%10s%18s%18s%24s%16s\n", "Size", "LCG_BaysDurham", "Xorshift128_MLT", "AdditiveLFG_RotateXor", "StandardMT");
    for (int sz : SIZES) {
        double tLCG, tXor, tLFG, tStd;
        {
            LCG_BaysDurham gen(baseSeed);
            tLCG = measureGenerationTime(gen, sz);
        }
        {
            Xorshift128_MLT gen(baseSeed + 1);
            tXor = measureGenerationTime(gen, sz);
        }
        {
            AdditiveLFG_RotateXor gen(baseSeed + 2);
            tLFG = measureGenerationTime(gen, sz);
        }
        {
            StandardMT gen(baseSeed + 3);
            tStd = measureGenerationTime(gen, sz);
        }

        printf("%10d%18.3f%18.3f%24.3f%16.3f\n", sz, tLCG, tXor, tLFG, tStd);
        fprintf(speed_f, "%d,%.3f,%.3f,%.3f,%.3f\n", sz, tLCG, tXor, tLFG, tStd);
    }
    fclose(speed_f);
    printf("\nSpeed results saved to output/speed_results.csv\n");

    return 0;
}