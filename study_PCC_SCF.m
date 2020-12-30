% Study this problem: is SC-Flip effective in decoding PCC-Polar
% codes? We should first approach the answer by Monte-Carlo simulations on
% the error distribution within one PC equation. If a partially decoded
% sequence is obtained, which ends with an erroneous parity bit, is the
% uncorrectly decoded bit identical to the bit with the most "uncertain"
% LLR?

% Three possibilities: identical, non-identical, >=3 (odd number) errors.
% construct SCFlip decoder.