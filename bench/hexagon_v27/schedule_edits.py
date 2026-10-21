q=O/'src/shared/affine_seam_discovery.hpp';s=q.read_text()
s=s.replace('size_t physical_round=0;double last_physical_time=elapsed();','size_t physical_round=0,last_physical_work=reducer.input_nonzeros();double last_physical_time=elapsed();')
s=s.replace('        last_physical_time=elapsed();','        last_physical_time=elapsed();last_physical_work=reducer.input_nonzeros();')
s=s.replace('if(elapsed()-last_physical_time>=physical_interval_seconds){','if(workload::schedule()?reducer.input_nonzeros()-last_physical_work>=16*options.batch_nonzeros:elapsed()-last_physical_time>=physical_interval_seconds){')
s=s.replace('bool stopped=false;size_t last_batch=reducer.completed_batches(),last_rank=reducer.rank(),stagnant=0;', '''bool stopped=false;size_t last_batch=reducer.completed_batches(),last_rank=reducer.rank(),stagnant=0;
        uint64_t last_arithmetic=reducer.arithmetic_work(),last_input=reducer.input_nonzeros();size_t last_fill=reducer.retained_nonzeros();''')
s=s.replace('size_t gain=reducer.rank()-last_rank;', '''size_t gain=reducer.rank()-last_rank;
                    uint64_t cost=reducer.arithmetic_work()-last_arithmetic+reducer.input_nonzeros()-last_input;
                    size_t fill=reducer.retained_nonzeros();if(fill>last_fill)cost+=fill-last_fill;
                    reducer.competing_work.observe(gain,cost);''')
s=s.replace('last_batch=reducer.completed_batches();last_rank=reducer.rank();', '''last_batch=reducer.completed_batches();last_rank=reducer.rank();
                    last_arithmetic=reducer.arithmetic_work();last_input=reducer.input_nonzeros();last_fill=reducer.retained_nonzeros();''')
q.write_text(s)
q=O/'src/shared/adaptive_physical_rows.hpp';s=q.read_text()
s=s.replace('    WorkCursor work;','''    WorkCursor work;
    structural_replay::Trace* trace=nullptr;
    std::vector<workload::GainWork> costs=std::vector<workload::GainWork>(3);
    workload::CycleProgress cycle_progress{2};
    uint64_t callback_work=0;size_t previous=SIZE_MAX,streak=0;''')
s=s.replace('size_t before=reducer.rank(),nnz=0,','uint64_t arithmetic_before=reducer.arithmetic_work();size_t before=reducer.rank(),nnz=0,')
s=s.replace('        size_t gain=reducer.rank()-before;','''        size_t gain=reducer.rank()-before;
        uint64_t measured_work=uint64_t(equations.size())*columns+nnz+reducer.arithmetic_work()-arithmetic_before;
        size_t kind=std::string(condition)=="collinear_tail_exact"?0:std::string(condition)=="mrk"?1:2;
        costs[kind].observe(gain,measured_work);callback_work+=measured_work;''')
s=s.replace('Policy value={})','Policy value={},structural_replay::Trace* plan=nullptr)')
s=s.replace('policy(value),columns(irr_physics::invariant_columns(context))','policy(value),columns(irr_physics::invariant_columns(context)),trace(plan)')
point=s.index('        run_budgeted_physics(')
s=s[:point]+'''        if(workload::schedule()){
            callback_work=0;
            // Recalibrate the integer budget from completed work, not elapsed
            // time. Cost and gain depend only on the supplied rows/field.
            uint64_t budget=4*std::min(costs[0].work,costs[1].work);
            budget=std::max(budget,uint64_t(columns)*policy.samples_per_batch);
            size_t tasks=0;
            while(nullity(reducer)>policy.target_augmented_nullity&&
                (trace&&trace->replay?trace->peek_decision()!=3:callback_work<budget&&tasks<policy.collinear_batches_per_cycle+1)){
                bool ope_due=!work.ope_done&&policy.ope_target_cap&&
                    (work.completed_cycles>=policy.cycles_before_ope||(collinear_done&&mrk_done));
                std::vector<bool> enabled={!collinear_done,!mrk_done,ope_due};
                // An early constraint is useful only if it competes with the
                // elimination it delays. Probe unmeasured tasks, then compare
                // observed rank gain per coefficient update, including fill.
                if(!(trace&&trace->replay))for(size_t i=0;i<enabled.size();++i)if(enabled[i]&&!workload::competitive(costs[i],reducer.competing_work)){
                    enabled[i]=false;std::cout<<"WORK_PHYSICS_DEFER round="<<round<<" task="<<i<<" physical_gain="<<costs[i].gain<<" physical_work="<<costs[i].work
                        <<" seam_gain="<<reducer.competing_work.gain<<" seam_work="<<reducer.competing_work.work<<std::endl;}
                size_t selected=enabled[2]?2:workload::choose(costs,enabled,previous,streak);
                if(selected==SIZE_MAX)break;
                if(trace)selected=trace->decision(selected);
                structural_replay::guard(selected<3&&enabled[selected],"Replay physical task no longer eligible");
                std::cout<<"WORK_PHYSICS_SELECT round="<<round<<" task="<<selected<<" consumed="<<callback_work<<" budget="<<budget<<" predicted_work="<<costs[selected].work<<" last_gain="<<costs[selected].gain<<std::endl;
                if(selected==2){apply_ope();work.ope_done=true;}
                else{batch(reducer,coordinate,selected==1,round,selected?policy.mrk_samples_per_batch:policy.samples_per_batch);
                    cycle_progress.admit(selected,{collinear_done,mrk_done});work.completed_cycles=cycle_progress.completed;}
                streak=selected==previous?streak+1:1;previous=selected;++tasks;
            }
            if(trace)structural_replay::guard(trace->decision(3)==3,"Replay physical stop changed");
            std::cout<<"WORK_PHYSICS_PAUSE round="<<round<<" work="<<callback_work<<" tasks="<<tasks<<" nullity="<<nullity(reducer)<<std::endl;
            return;
        }
''' +s[point:]
q.write_text(s)
