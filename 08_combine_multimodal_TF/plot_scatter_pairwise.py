import seaborn as sns
import matplotlib.pyplot as plt
import scipy


def plot_scatter_signature_tf(df, x_column, y_column, outpath):
    sns.relplot(y=y_column,
                x=x_column,
                data=df,
                hue='Sample Type',
                ax=ax,
                palette=sns.color_palette(["black", "#DD6FA3", "#9D58F4", "#EC9D40", ]),
                )
    plt.xlim(-0.03, 1.1)
    plt.ylim(-0.03, 1.1)

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(df[x_column], df[y_column])
    # slope, intercept = np.polyfit(tf_df["ILL_TF_ichorCNA"],tf_df["CYC_TF_ichorCNA"], 1)

    # add linear regression line to scatterplot
    plt.plot([0, 1], [intercept, slope * 1 + intercept], color='grey')
    # math symbol between $$ signs. https://matplotlib.org/2.0.2/users/mathtext.html
    plt.text(0.7, 0.1, f"$R^2$={round(r_value ** 2, 4)}")
    print(f'y={round(slope, 3)}X+{round(intercept, 3)}', )
    print("r-square=", round(r_value ** 2, 4))
    plt.savefig(outpath, dpi=300)
    plt.savefig(outpath + '.png', dpi=300)
    plt.show()