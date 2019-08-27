// void fullHaploTest() {
//     Haplotyper h(5, 2, 33, 8, 5, 3, 5, false, true, Options::FilterType::GAUSS);
//
//     equals(h.getOffset(), 10);
//
//     // First 5 pushes inside basespreader
//     for (int i = 0; i < 5; ++i) {
//         equals(h.push("C", "}", 0, 'A'), 0);
//     }
//
//     // next 11 pushes inside filterbuffer approaching 0.73
//     for (int i = 0; i < 10; ++i) {
//         h.push("C", "}", 0, 'A');
//     }
//     equals(h.push("C", "}", 0, 'A'), 7);
//     equals(h.push("C", "}", 0, 'A'), 7);
//
//     // reset
//     Haplotyper h2(5, 2, 33, 8, 5, 3, 5, false, true, Options::FilterType::GAUSS);
//
//     h2.push("CCC", "}}}", 15, 'A');
//
//     // push spreaded bases into filterbuffer
//     for (int i = 0; i < 8; ++i) {
//         h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C');  // Activity close to 0
//     }
//     // Reach maximum
//     equals(h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C'), 7);
//     equals(h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C'), 7);
//     equals(h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C'), 7);
//
//     // Reach minimum
//     for (int i = 0; i < 9; ++i) {
//         h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C');
//     }
//
//     equals(h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C'), 0);
// }
