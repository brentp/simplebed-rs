use criterion::{black_box, criterion_group, criterion_main, Criterion};
use simplebed::BedRecord;

pub fn bench_parsing(c: &mut Criterion) {
    let mut group = c.benchmark_group("parse_line");

    group.bench_function("simple_line", |b| {
        let line = "chr1\t1000\t2000";
        b.iter(|| BedRecord::parse_line(black_box(line)).unwrap())
    });

    group.bench_function("full_line", |b| {
        let line = "chr1\t1000\t2000\tgene1\t100.5\tstring_field\t42\t3.14";
        b.iter(|| BedRecord::parse_line(black_box(line)).unwrap())
    });

    group.bench_function("long_chrom", |b| {
        let line = "chr1_very_long_contig_name_that_might_slow_things_down\t1000\t2000\tgene1";
        b.iter(|| BedRecord::parse_line(black_box(line)).unwrap())
    });

    group.finish();
}

criterion_group!(parse_line_benches, bench_parsing);
criterion_main!(parse_line_benches);
