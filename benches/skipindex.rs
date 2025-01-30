use criterion::{black_box, criterion_group, criterion_main, Criterion};
use simplebed::{BedReader, BedRecord};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

fn setup_test_bed() -> std::io::Result<()> {
    let test_dir = Path::new("benches/data");
    if !test_dir.exists() {
        std::fs::create_dir_all(test_dir)?;
    }

    let path = test_dir.join("bench.bed");
    let mut file = BufWriter::new(File::create(path)?);

    // Write entries for each chromosome from chr1 to chr22
    for chrom in 1..=22 {
        let chrom_name = format!("chr{}", chrom);

        for i in 0..50_000 {
            // 10x more entries than test for better benchmarking
            let start = i * 100;
            let end = start + 50;
            let score = (i % 100) as f64;
            let name = format!("feat_{}_{}", chrom, i);

            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}",
                chrom_name, start, end, name, score,
            )?;
        }
    }
    file.flush()?;
    Ok(())
}

fn get_chromosome_order() -> HashMap<String, usize> {
    (1..=22)
        .map(|i| (format!("chr{}", i), i - 1))
        .collect::<HashMap<_, _>>()
}

fn bench_index_creation(c: &mut Criterion) {
    setup_test_bed().expect("Failed to create test bed file");
    let path = Path::new("benches/data/bench.bed");

    c.bench_function("create_skip_index", |b| {
        b.iter(|| {
            let mut bed_reader = BedReader::<File>::from_path(path).unwrap();
            bed_reader.set_chromosome_order(get_chromosome_order());
            bed_reader.build_skip_index().unwrap();
        })
    });
}

fn bench_queries(c: &mut Criterion) {
    setup_test_bed().expect("Failed to create test bed file");
    let path = Path::new("benches/data/bench.bed");

    // Setup reader with index
    let mut bed_reader = BedReader::<File>::from_path(path).expect("Failed to open bed file");
    bed_reader.set_chromosome_order(get_chromosome_order());
    bed_reader
        .build_skip_index()
        .expect("Failed to build index");

    let mut group = c.benchmark_group("queries");

    // Benchmark small region query
    group.bench_function("small_region_query", |b| {
        b.iter(|| {
            let records: Vec<BedRecord> = bed_reader
                .query("chr10", 11000, 11200)
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            black_box(records)
        })
    });

    // Benchmark chromosome switch
    group.bench_function("chromosome_switch_query", |b| {
        b.iter(|| {
            for chrom in ["chr5", "chr15", "chr7", "chr20"] {
                let records: Vec<BedRecord> = bed_reader
                    .query(chrom, 11000, 11200)
                    .unwrap()
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap();
                black_box(records);
            }
        })
    });

    group.finish();
}

criterion_group!(benches, bench_index_creation, bench_queries);
criterion_main!(benches);
