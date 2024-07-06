module sv

using FileIO
using Plots
using Gtk.ShortNames,Gtk

export save_variables, save_plots, select_folder

# Function to save variables
function save_variables(filename, vars...)
    open(filename, "w") do file
        varname = ["kx", "δ", "A", "cloud_r", "cloud_T","n"]
        for (i, var) in enumerate(vars)
            if i<=length(varname)
                println(file, "$(varname[i]): $var")
            else
                println(file,"$var")
            end
        end
    end
end

# Function to save plots
function save_plots(folder, plots...)
    name = ["wavefront","2d","1d_x","1d_y"]
    for (i, plot) in enumerate(plots)
        savefig(plot, joinpath(folder, "$(name[i]).png"))
    end
end

# Function to select folder
function select_folder()
    dialog = GtkFileChooserDialog("Select Folder", nothing, GtkFileChooserAction.SELECT_FOLDER,
                                  "Cancel", GtkResponseType(:CANCEL),
                                  "Select", GtkResponseType(:ACCEPT))
    if GtkDialog.run(dialog) == GtkResponseType(:ACCEPT)
        folder = GtkFileChooser.get_filename(dialog)
        GtkWidget.destroy(dialog)
        return folder
    else
        GtkWidget.destroy(dialog)
        return nothing
    end
end


end # module SaveUtils
