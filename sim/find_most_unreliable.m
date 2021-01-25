function ret_index = find_most_unreliable(is_used, LLRs)
    LLRs_unused = LLRs(~is_used);
    [~,ind] = min(LLRs_unused);
    tmp = find(~is_used, ind);
    ret_index = tmp(end);
end